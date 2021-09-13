#include "timeserieslistport.h"

namespace voreen {
    TimeSeriesListPort::TimeSeriesListPort(PortDirection direction, const std::string& id, const std::string& guiName,
        bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
        : GenericPort<TimeSeriesList>(direction, id, guiName, allowMultipleConnections, invalidationLevel) {

    }

    std::string TimeSeriesListPort::getClassName() const {
        return "TimeSeriesListPort";
    }
    Port* TimeSeriesListPort::create(PortDirection direction, const std::string& id, const std::string& guiName) const {
        return new TimeSeriesListPort(direction, id, guiName);
    }
    tgt::col3 TimeSeriesListPort::getColorHint() const {
        return tgt::col3(135.15, 40.8, 81.6); // plum pudding
    }
    std::string TimeSeriesListPort::getContentDescription() const {
        std::stringstream strstr;
        strstr << Port::getContentDescription();
        if (hasData()) {
            const voreen::TimeSeriesList* timeSeriesList = getData();
            strstr << std::endl << "Number of time series: " << std::to_string(timeSeriesList->getNumberOfTimeSeries());
            strstr << std::endl << "Overall start: " << std::to_string(timeSeriesList->getOverallStartTime());
            strstr << std::endl << "Overall end: " << std::to_string(timeSeriesList->getOverallEndTime());
            
            std::string fields; 
            if (timeSeriesList->getNumberOfFields() > 0) 
                fields = timeSeriesList->getFieldName(0);
            for (size_t i = 1; i < timeSeriesList->getNumberOfFields(); ++i) {
                fields += ", " + timeSeriesList->getFieldName(i);
            }
            strstr << std::endl << "Field names: " << fields;

            std::string components;
            if (timeSeriesList->getNumberOfComponentsPerField() > 0)
                components = std::to_string(timeSeriesList->getComponentForField(0));
            for (size_t i = 1; i < timeSeriesList->getNumberOfComponentsPerField(); ++i) {
                components += ", " + std::to_string(timeSeriesList->getComponentForField(i));
            }
            strstr << std::endl << "Components per field: " << components;
        }
        return strstr.str();
    }
    std::string TimeSeriesListPort::getContentDescriptionHTML() const {
        std::stringstream strstr;
        strstr << Port::getContentDescriptionHTML();
        if (hasData()) {
            const voreen::TimeSeriesList* timeSeriesList = getData();
            strstr << "<br>" << "Number of time series: " << std::to_string(timeSeriesList->getNumberOfTimeSeries());
            strstr << "<br>" << "Overall start: " << std::to_string(timeSeriesList->getOverallStartTime());
            strstr << "<br>" << "Overall end: " << std::to_string(timeSeriesList->getOverallEndTime());
            
            std::string fields;
            if (timeSeriesList->getNumberOfFields() > 0)
                fields = timeSeriesList->getFieldName(0);
            for (size_t i = 1; i < timeSeriesList->getNumberOfFields(); ++i) {
                fields += ", " + timeSeriesList->getFieldName(i);
            }
            strstr << "<br>" << "Field names: " << fields;

            std::string components;
            if (timeSeriesList->getNumberOfComponentsPerField() > 0)
                components = std::to_string(timeSeriesList->getComponentForField(0));
            for (size_t i = 1; i < timeSeriesList->getNumberOfComponentsPerField(); ++i) {
                components += ", " + std::to_string(timeSeriesList->getComponentForField(i));
            }
            strstr << "<br>" << "Components per field: " << components;
        }
        return strstr.str();
    }
}