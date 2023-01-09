import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import plotly.express as px

import voreen
import json
import time
import os
from threading import Thread
from multiprocessing import Process, Manager

# We create a global proxy list that can be copied and used from within any process and scope.
# It is used to communicate between processes.
manager = Manager()
from_voreen_updates_lock = manager.Lock()
from_voreen_updates = manager.dict()
from_dash_updates_lock = manager.Lock()
from_dash_updates = manager.dict()

def get_property_key(processor_id, property_id):
    return f'{processor_id}: {property_id}'

# The Voreen Dash App
class VoreenDashApp(Thread):
    def __init__(self, target):
        super().__init__()

        self.server = Process(target=target)
        self.server_terminated = False
        self.start_server()

        self.setDaemon(True)
        self.start()

        while not self.server_terminated:
            time.sleep(1)
            pass

    def start_server(self):
        self.server.start()

    def terminate_server(self):
        self.server.terminate()
        self.server.join()
        self.server_terminated = True

    def run(self):
        while True:

            if voreen.exitRequested():
                self.terminate_server()
                break

            time.sleep(1)


class VoreenPropertyWatcher:
    def __init__(self, processor_id, property_id):
        self.processor_id = processor_id
        self.property_id = property_id
        self.property_value = voreen.getPropertyValue(processor_id, property_id)
        self.property_value_min = voreen.getPropertyMinValue(processor_id, property_id)
        self.property_value_max = voreen.getPropertyMaxValue(processor_id, property_id)

    def serialize(self, value):
        return value  # TODO: serialization might be necessary depending on datatype
        try:
            s = json.dumps(value.__dict__)
        except Exception as e:
            s = json.dump(value)
        return s

    def compare(self):
        return self.serialize(voreen.getPropertyValue(self.processor_id, self.property_id)) == self.serialize(self.property_value)

    def fetch(self):
        if not self.compare():
            self.property_value = voreen.getPropertyValue(self.processor_id, self.property_id)
            key = get_property_key(self.processor_id, self.property_id)
            from_voreen_updates_lock.acquire()
            from_voreen_updates[key] = self.property_value
            from_voreen_updates_lock.release()

            return True

        return False


class VoreenPropertyWatcherThread(Thread):
    def __init__(self):
        super().__init__()
        self.watchers = {}
        self.setDaemon(True)
        self.start()

    def add_watch(self, processor_id, property_id):
        key = get_property_key(processor_id, property_id)
        self.watchers[key] = VoreenPropertyWatcher(processor_id, property_id)

    def remove_watch(self, processor_id, property_id):
        key = get_property_key(processor_id, property_id)
        self.watchers.pop(key)

    def run(self):
        while (True):

            for k in self.watchers:
                self.watchers[k].fetch()

            from_dash_updates_lock.acquire()
            for k in from_dash_updates:
                self.watchers[k].property_value = from_dash_updates[k][-1]
                voreen.setPropertyValue(*from_dash_updates[k])
            from_dash_updates.clear()
            from_dash_updates_lock.release()

            time.sleep(0.1)


def setup_dash_app():
    app = dash.Dash()

    def create_slider_widget(watch):
        key = get_property_key(watch.processor_id, watch.property_id)
        div_key = key+'_container'
        interval_key = key+'_interval'
        dummy = html.Div(id=div_key)
        label = html.H5(key)
        slider = dcc.Slider(id=key,
                            min=watch.property_value_min,
                            max=watch.property_value_max,
                            value=watch.property_value,
#                            step=(watch.property_value_max-watch.property_value_min)/50,
                            updatemode='drag')
        interval = dcc.Interval(id=interval_key, interval=100, n_intervals=0)

        @app.callback(
            Output(div_key, 'children'),
            Input(key, 'value')
        )
        def update_slider(value):
            from_dash_updates_lock.acquire()
            from_dash_updates[key] = (watch.processor_id, watch.property_id, value)
            from_dash_updates_lock.release()
            return None

        @app.callback(
            Output(key, 'value'),
            Input(interval_key, 'n_intervals'),
            State(key, 'value')
        )
        def update_interval(_, value):
            from_voreen_updates_lock.acquire()
            if key in from_voreen_updates:
                value = from_voreen_updates[key]
                from_voreen_updates.pop(key)
                from_voreen_updates_lock.release()
                return value
            from_voreen_updates_lock.release()
            raise dash.exceptions.PreventUpdate

        return dummy, key, slider, interval

    def create_widgets_from_properties():
        divs = []
        for key in watcher.watchers:
            divs.extend(create_slider_widget(watcher.watchers[key]))

        return divs

    app.layout = html.Div(children=create_widgets_from_properties())
    app.run_server(debug=False, port=8050, use_reloader=False)


# Setup watcher!
watcher = VoreenPropertyWatcherThread()
watcher.add_watch('IsosurfaceExtractor', 'isoValueProp')
watcher.add_watch('Background', 'angle')

# Setup Voreen Dash App!
VoreenDashApp(setup_dash_app)

