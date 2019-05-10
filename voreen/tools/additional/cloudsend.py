#!/usr/bin/env python3
import requests
import os
import argparse
from xml.etree import ElementTree
import threading
from collections import namedtuple
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--server', help='Server URL', type=str, default='https://uni-muenster.sciebo.de')
parser.add_argument('-t', '--share-folder-token', help='Share folder token on the server', type=str, default='YdhJ4czrPsfqPL6')
parser.add_argument('-p', '--password', help='filter a parameter (<parametername>=<value>)', type=str)
parser.add_argument('-u', '--upload-dir', help='upload the specified directory recursively to the destination', type=str)
parser.add_argument('--max-old-runs', help='maximum number of runs to keep when adding a ci report', type=int, default=3)
parser.add_argument('-d', '--delete', help='delete destination', action="store_true")
parser.add_argument('-f', '--file', help='upload file as destionation', type=str)
parser.add_argument('-l', '--list', help='list remote directory', action="store_true")
parser.add_argument('--add-ci-report', help='Add the provided ci report directory')
parser.add_argument('destination')
args = parser.parse_args()

pubsuffix='public.php/webdav'
headers = {'X-Requested-With': 'XMLHttpRequest'}
reports_web_server_root = "https://sso.uni-muenster.de/IVV5WS/vrreports"
overview_file_name = "overview.html"
overview_css_name = "overview.css"
regressiontest_report_file = "regressiontest-report.html"

url_base = args.server + '/' + pubsuffix

auth = (args.share_folder_token, args.password)

max_server_error_retries=5

http_request_timeout_seconds=10

def retry_request(request):
    code = 500
    num_tries = 0
    rep = request()
    while rep.status_code == 500 and num_tries < max_server_error_retries:
        try:
            rep = request()
        except requests.exceptions.ReadTimeout as e:
            print("Request timeout: {}", e)
        num_tries += 1
    return rep

def put_as_file(data, destination):
    url = url_base + '/' + destination
    r = retry_request(lambda: requests.put(url, data=data, headers=headers, auth=auth, timeout=http_request_timeout_seconds))

    if r.status_code < 200 or r.status_code >= 300:
        print('Sending file {} not successful: {} {}'.format(destination, r.status_code, r.text))
    #print('Sending file {}\nReply: {} {}'.format(destination, r.status_code, r.text))
def put_single_file(file_to_send, destination):
    with open(file_to_send, 'rb') as f:
        put_as_file(f.read(), destination)

def make_dir(destination):
    url = url_base + '/' + destination
    r = retry_request(lambda: requests.request('MKCOL', url, headers=headers, auth=auth, timeout=http_request_timeout_seconds))
    if r.status_code < 200 or r.status_code >= 300:
        print('Creating dir {} not successful: {} {}'.format(destination, r.status_code, r.text))

def ensure_path_exists(destination):
    existing_path = ""
    for part in destination.split('/'):
        make_dir(existing_path + '/' + part)
        existing_path += '/' + part

def put_tree(dir_path, destination):
    ensure_path_exists(destination)
    root_len = len(dir_path)
    for root, dirs, files in os.walk(dir_path, topdown=True):
        threads = []
        root_relative_path = root[root_len:]
        destination_dir_prefix = destination + '/' + root_relative_path + '/'
        for d in dirs:
            t = threading.Thread(target=make_dir, args=(destination_dir_prefix +  d,))
            #make_dir(destination_dir_prefix +  d)
            t.start()
            threads.append(t)
        for f in files:
            t = threading.Thread(target=put_single_file, args=(os.path.join(root, f), destination_dir_prefix + f))
            t.start()
            threads.append(t)
            #put_single_file(os.path.join(root, f), destination_dir_prefix + f)

        for t in threads:
            t.join()

def delete(destination):
    url = url_base + '/' + destination
    r = retry_request(lambda: requests.delete(url, headers=headers, auth=auth, timeout=http_request_timeout_seconds))
    if r.status_code < 200 or r.status_code >= 300:
        print('Deleting {} not successful: {} {}'.format(destination, r.status_code, r.text))

def list_dir(destination):
    data = '''<?xml version="1.0" encoding="UTF-8"?>
     <d:propfind xmlns:d="DAV:">
       <d:prop xmlns:oc="http://owncloud.org/ns">
         <d:resourcetype/>
       </d:prop>
     </d:propfind>'''
    url = url_base + '/' + destination
    r = retry_request(lambda: requests.request('PROPFIND', url, data=data, headers=headers, auth=auth, timeout=http_request_timeout_seconds))
    if r.status_code < 200 or r.status_code >= 300:
        print('Listing dir {} not successful: {} {}'.format(destination, r.status_code, r.text))
        sys.exit(1)

    tree = ElementTree.fromstring(r.content)
    shortest_len = None
    candidates = []
    for node in tree.findall('ns0:response/ns0:href', namespaces={
        "ns0": "DAV:"
        }):
        text = node.text
        if shortest_len is None or shortest_len > len(text):
            shortest_len = len(text)
        candidates.append(text)
    listing = []
    for candidate in candidates:
        shortened = candidate[shortest_len:]
        if shortened:
            listing.append(shortened)

    return listing


class Run(object):
    def path(self):
        return "{}--{}--{}--{}".format(self.branch, self.date, self.commit, self.build)
    def parse(line):
        components = line.split("--")
        if len(components) < 4:
            return None
        # Allow for "--" in branch name, but not in build (or commit hash or time)
        run = Run()
        run.build = components[-1]
        run.commit = components[-2]
        run.date = components[-3]
        run.branch = "--".join(components[:-3])
        return run
    def __repr__(self):
        return self.path()
    def __str__(self):
        return self.path()

if args.list:
    print(list_dir(args.destination))
if args.delete:
    delete(args.destination)
if args.upload_dir:
    put_tree(args.upload_dir, args.destination)
if args.file:
    put_single_file(args.file, args.destination)
if args.add_ci_report:
    report_base_dir = os.path.dirname(args.destination)
    report_name = os.path.basename(args.destination)
    run = Run.parse(report_name)
    if run is None:
        print("Invalid run format.\nExpected: '<branch>--<datetime>--<commit>--<build>'\nGot: '{}'".format(args.destination))
        sys.exit(1)

    print("Adding run: branch: {}, time: {}, commit: {}, build: {}".format(run.branch, run.date, run.commit, run.build))
    put_tree(args.add_ci_report, args.destination)

    listing = list_dir("/" + report_base_dir)
    runs = []
    matching_runs = []
    for s in listing:
        r = Run.parse(os.path.normpath(s))
        if r:
            runs.append(r)
            if r.branch == run.branch and r.build == run.build:
                matching_runs.append(r)

    matching_runs.sort(key=lambda r: r.date, reverse=True)
    to_delete = matching_runs[args.max_old_runs:]
    for d in to_delete:
        print("Deleting old run: {}".format(d.path()))
        delete(report_base_dir + "/" + d.path())
        runs.remove(d)

    runs.sort(key=lambda r: r.date, reverse=True)
    index_page  = '<!DOCTYPE html>'
    index_page += '<html>'
    index_page += '<head><link rel="stylesheet" href="{}"></head>\n'.format(overview_css_name)
    index_page += '<body>\n'
    index_page += '<marquee>Voreen Regression Test Reports</marquee>\n'
    branches = list(set([r.branch for r in runs]))
    branches.sort()
    for branch in branches:
        index_page += '<table><tbody>\n'
        index_page += '<tr><td class="branchheader" colspan="3"><b>{}</b></td></tr>\n'.format(branch)
        index_page += '<th>Commit</th><th>Time</th><th>Build</th>\n'
        for r in runs:
            if r.branch == branch:
                index_page += '<tr><td class="commit"><a href="{}/{}/{}">{}</a></td><td class="date">{}</td><td class="build">{}</td></tr>\n'.format(report_base_dir, r.path(), regressiontest_report_file, r.commit, r.date, r.build)
        index_page += '</tbody></table>\n'
    index_page += '</body>\n'
    index_page += '</html>'

    put_as_file(index_page, overview_file_name)
    print("")
    print("New testreport: {}/{}/{}".format(reports_web_server_root, run.path(), regressiontest_report_file))
    print("Overview: {}/{}".format(reports_web_server_root, overview_file_name))
