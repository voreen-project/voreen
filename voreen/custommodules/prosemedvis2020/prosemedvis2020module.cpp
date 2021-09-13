#include "prosemedvis2020module.h"

// MRI
#include "processors/concretevesselgraphcreator.h"
#include "processors/labelvolumecreator.h"
#include "processors/templatevesselgraphloader.h"
#include "processors/concretevesselgraphsource.h"
#include "processors/abstractvisualizationcreator.h"
#include "processors/cvglabelswapper.h"

// PET
#include "processors/timeseriesextraction.h"
#include "processors/meanplot.h"
#include "processors/ensemblevolumepicking.h"
#include "processors/clickabletextureoverlay.h"
#include "processors/timeseriesfilter.h"

namespace voreen {

ProseMedvis2020Module::ProseMedvis2020Module(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("ProSeMedVis2020Module");

    setGuiName("ProSe MedVis 2020 Module");

    // MRI
    registerProcessor(new ConcreteVesselGraphCreator());
    registerProcessor(new LabelVolumeCreator());
    registerProcessor(new TemplateVesselGraphLoader());
    registerProcessor(new ConcreteVesselGraphSource());
    registerProcessor(new AbstractVisualizationCreator());
    registerProcessor(new CvgLabelSwapper());

    // PET
    registerProcessor(new TimeseriesExtraction());
    registerProcessor(new MeanPlot());
    registerProcessor(new EnsembleVolumePicking());
    registerProcessor(new ClickableTextureOverlay());
    registerProcessor(new TimeseriesFilter());

    addShaderPath(getModulePath("glsl"));
}

}
