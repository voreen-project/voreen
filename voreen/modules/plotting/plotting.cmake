# external dependency: triangle library
SET(MOD_INCLUDE_DIRECTORIES
    ${MOD_DIR}/ext/triangle
)
SET(MOD_INSTALL_FILES
    ${MOD_DIR}/ext/triangle/license.txt
)

# sample data
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/data
)

################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS PlottingModule)

SET(MOD_CORE_SOURCES
    # processor sources
    ${MOD_DIR}/processors/barplot.cpp
    ${MOD_DIR}/processors/hemisphereplot.cpp
    ${MOD_DIR}/processors/imageanalyzer.cpp
    ${MOD_DIR}/processors/lineplot.cpp
    ${MOD_DIR}/processors/plotdatafitfunction.cpp
    ${MOD_DIR}/processors/plotdatagroup.cpp
    ${MOD_DIR}/processors/plotdatasource.cpp
    ${MOD_DIR}/processors/plotdataexport.cpp
    ${MOD_DIR}/processors/plotdataexporttext.cpp
    ${MOD_DIR}/processors/plotdataexportbase.cpp
    ${MOD_DIR}/processors/plotdataselect.cpp
    ${MOD_DIR}/processors/plotdatamerge.cpp
    ${MOD_DIR}/processors/plotfunctiondiscret.cpp
    ${MOD_DIR}/processors/plotfunctionsource.cpp
    ${MOD_DIR}/processors/plotprocessor.cpp
    ${MOD_DIR}/processors/surfaceplot.cpp
    ${MOD_DIR}/processors/scatterplot.cpp

    # Property sources
    ${MOD_DIR}/properties/colormapproperty.cpp
    ${MOD_DIR}/properties/plotentitiesproperty.cpp
    ${MOD_DIR}/properties/plotpredicateproperty.cpp
    ${MOD_DIR}/properties/plotselectionproperty.cpp
    ${MOD_DIR}/properties/plotdataproperty.cpp
    ${MOD_DIR}/properties/link/linkevaluatorplotselection.cpp

    # Port sources
    ${MOD_DIR}/ports/plotport.cpp

    # Datastructures sources
    ${MOD_DIR}/datastructures/aggregationfunction.cpp
    ${MOD_DIR}/datastructures/colormap.cpp
    ${MOD_DIR}/datastructures/plotbase.cpp
    ${MOD_DIR}/datastructures/plotcell.cpp
    ${MOD_DIR}/datastructures/plotdata.cpp
    ${MOD_DIR}/datastructures/plotentitysettings.cpp
    ${MOD_DIR}/datastructures/plotexpression.cpp
    ${MOD_DIR}/datastructures/plotfunction.cpp
    ${MOD_DIR}/datastructures/plotpredicate.cpp
    ${MOD_DIR}/datastructures/plotrow.cpp
    ${MOD_DIR}/datastructures/plotselection.cpp

    # Utils sources
    ${MOD_DIR}/utils/functionlibrary.cpp
    ${MOD_DIR}/utils/plotdatainserter.cpp
    ${MOD_DIR}/utils/plotlibrary/plotlibrary.cpp
    ${MOD_DIR}/utils/plotlibrary/plotlibraryfilebase.cpp
    ${MOD_DIR}/utils/plotlibrary/plotlibrarylatex.cpp
    ${MOD_DIR}/utils/plotlibrary/plotlibrarylatexrender.cpp
    ${MOD_DIR}/utils/plotlibrary/plotlibrarynonefilebase.cpp
    ${MOD_DIR}/utils/plotlibrary/plotlibraryopengl.cpp
    ${MOD_DIR}/utils/plotlibrary/plotlibrarysvg.cpp
    ${MOD_DIR}/utils/plotlibrary/plotlibrarysvgrender.cpp
    ${MOD_DIR}/utils/plotpredicatefactory.cpp
    ${MOD_DIR}/utils/smartlabel.cpp
    ${MOD_DIR}/utils/parser/plotfunctionterminals.cpp
    ${MOD_DIR}/utils/parser/plotfunctiongrammar.cpp
    ${MOD_DIR}/utils/parser/plotfunctionlexer.cpp
    ${MOD_DIR}/utils/parser/plotfunctionparser.cpp
    ${MOD_DIR}/utils/parser/plotfunctionnode.cpp
    ${MOD_DIR}/utils/parser/plotfunctionvisitor.cpp
    ${MOD_DIR}/utils/parser/plotfunctiontoken.cpp

    # Interaction sources
    ${MOD_DIR}/interaction/plotcamerainteractionhandler.cpp
    ${MOD_DIR}/interaction/plotpickingmanager.cpp

    # ext sources
    ${MOD_DIR}/ext/triangle/src/del_impl.cpp
)

SET(MOD_CORE_HEADERS
    # processor headers
    ${MOD_DIR}/processors/plotdatasource.h
    ${MOD_DIR}/processors/plotfunctionsource.h
    ${MOD_DIR}/processors/barplot.h
    ${MOD_DIR}/processors/hemisphereplot.h
    ${MOD_DIR}/processors/imageanalyzer.h
    ${MOD_DIR}/processors/lineplot.h
    ${MOD_DIR}/processors/plotdatafitfunction.h
    ${MOD_DIR}/processors/plotdatagroup.h
    ${MOD_DIR}/processors/plotdataexport.h
    ${MOD_DIR}/processors/plotdataexporttext.h
    ${MOD_DIR}/processors/plotdataexportbase.h
    ${MOD_DIR}/processors/plotdataselect.h
    ${MOD_DIR}/processors/plotdatamerge.h
    ${MOD_DIR}/processors/plotfunctiondiscret.h
    ${MOD_DIR}/processors/plotprocessor.h
    ${MOD_DIR}/processors/surfaceplot.h
    ${MOD_DIR}/processors/scatterplot.h

    # Property headers
    ${MOD_DIR}/properties/colormapproperty.h
    ${MOD_DIR}/properties/plotentitiesproperty.h
    ${MOD_DIR}/properties/plotpredicateproperty.h
    ${MOD_DIR}/properties/plotselectionproperty.h
    ${MOD_DIR}/properties/plotdataproperty.h
    ${MOD_DIR}/properties/link/linkevaluatorcolormapid.h
    ${MOD_DIR}/properties/link/linkevaluatorplotselection.h
    ${MOD_DIR}/properties/link/linkevaluatorplotentitiesid.h

    # Port headers
    ${MOD_DIR}/ports/plotport.h

    # Datastructures headers
    ${MOD_DIR}/datastructures/aggregationfunction.h
    ${MOD_DIR}/datastructures/colormap.h
    ${MOD_DIR}/datastructures/interval.h
    ${MOD_DIR}/datastructures/plotbase.h
    ${MOD_DIR}/datastructures/plotcell.h
    ${MOD_DIR}/datastructures/plotdata.h
    ${MOD_DIR}/datastructures/plotentitysettings.h
    ${MOD_DIR}/datastructures/plotexpression.h
    ${MOD_DIR}/datastructures/plotfunction.h
    ${MOD_DIR}/datastructures/plotpredicate.h
    ${MOD_DIR}/datastructures/plotrow.h
    ${MOD_DIR}/datastructures/plotselection.h
    ${MOD_DIR}/datastructures/plotzoomstate.h

    # Utils headers
    ${MOD_DIR}/utils/functionlibrary.h
    ${MOD_DIR}/utils/plotdatainserter.h
    ${MOD_DIR}/utils/plotlibrary/plotlibrarylatex.h
    ${MOD_DIR}/utils/plotlibrary/plotlibrarylatexrender.h
    ${MOD_DIR}/utils/plotlibrary/plotlibrary.h
    ${MOD_DIR}/utils/plotlibrary/plotlibraryfilebase.h
    ${MOD_DIR}/utils/plotlibrary/plotlibrarynonefilebase.h
    ${MOD_DIR}/utils/plotlibrary/plotlibraryopengl.h
    ${MOD_DIR}/utils/plotlibrary/plotlibrarysvg.h
    ${MOD_DIR}/utils/plotlibrary/plotlibrarysvgrender.h
    ${MOD_DIR}/utils/plotpredicatefactory.h
    ${MOD_DIR}/utils/smartlabel.h
    ${MOD_DIR}/utils/parser/plotfunctionterminals.h
    ${MOD_DIR}/utils/parser/plotfunctiongrammar.h
    ${MOD_DIR}/utils/parser/plotfunctionlexer.h
    ${MOD_DIR}/utils/parser/plotfunctionparser.h
    ${MOD_DIR}/utils/parser/plotfunctionnode.h
    ${MOD_DIR}/utils/parser/plotfunctionvisitor.h
    ${MOD_DIR}/utils/parser/plotfunctiontoken.h

    # Interaction headers
    ${MOD_DIR}/interaction/plotcamerainteractionhandler.h
    ${MOD_DIR}/interaction/plotpickingmanager.h

    # ext headers
    ${MOD_DIR}/ext/triangle/include/del_interface.hpp
    ${MOD_DIR}/ext/triangle/include/dpoint.hpp
    ${MOD_DIR}/ext/triangle/include/triangle.h
    ${MOD_DIR}/ext/triangle/include/triangle_impl.hpp
)


################################################################################
# Qt module resources
################################################################################
SET(MOD_QT_MODULECLASS PlottingModuleQt)

SET(MOD_QT_SOURCES
    # widget sources
    ${MOD_DIR}/qt/plotdatawidget.cpp
    ${MOD_DIR}/qt/plotdatasimpletablemodel.cpp
    ${MOD_DIR}/qt/extendedtable.cpp
    ${MOD_DIR}/qt/plotdataextendedtablemodel.cpp
    ${MOD_DIR}/qt/plotentitysettingsdialog.cpp
    ${MOD_DIR}/qt/plotpredicatedialog.cpp
    ${MOD_DIR}/qt/plotselectionentrytablemodel.cpp
    ${MOD_DIR}/qt/plotselectiontablemodel.cpp
    ${MOD_DIR}/qt/plotselectiondialog.cpp

    # Processor widget sources
    ${MOD_DIR}/qt/processorwidgets/plotprocessorwidgetfactory.cpp
    ${MOD_DIR}/qt/processorwidgets/plotdatamergewidget.cpp
    ${MOD_DIR}/qt/processorwidgets/plotdatafitfunctionwidget.cpp
    ${MOD_DIR}/qt/processorwidgets/plotdataselectwidget.cpp
    ${MOD_DIR}/qt/processorwidgets/plotdatagroupwidget.cpp
    ${MOD_DIR}/qt/processorwidgets/plotdataprocessorwidget.cpp

    # Property widget sources
    ${MOD_DIR}/qt/propertywidgets/colormappropertywidget.cpp
    ${MOD_DIR}/qt/propertywidgets/plotpropertywidgetfactory.cpp
    ${MOD_DIR}/qt/propertywidgets/plotdatapropertywidget.cpp
    ${MOD_DIR}/qt/propertywidgets/plotentitiespropertywidget.cpp
    ${MOD_DIR}/qt/propertywidgets/plotpredicatepropertywidget.cpp
    ${MOD_DIR}/qt/propertywidgets/plotselectionpropertywidget.cpp
)

SET(MOD_QT_HEADERS
    # Widget headers
    ${MOD_DIR}/qt/plotdatawidget.h
    ${MOD_DIR}/qt/plotdatasimpletablemodel.h
    ${MOD_DIR}/qt/extendedtable.h
    ${MOD_DIR}/qt/plotdataextendedtablemodel.h
    ${MOD_DIR}/qt/plotentitysettingsdialog.h
    ${MOD_DIR}/qt/plotpredicatedialog.h
    ${MOD_DIR}/qt/plotselectionentrytablemodel.h
    ${MOD_DIR}/qt/plotselectiontablemodel.h
    ${MOD_DIR}/qt/plotselectiondialog.h

    # Processor widget headers
    ${MOD_DIR}/qt/processorwidgets/plotdatamergewidget.h
    ${MOD_DIR}/qt/processorwidgets/plotdatafitfunctionwidget.h
    ${MOD_DIR}/qt/processorwidgets/plotdataselectwidget.h
    ${MOD_DIR}/qt/processorwidgets/plotdatagroupwidget.h
    ${MOD_DIR}/qt/processorwidgets/plotdataprocessorwidget.h

    # Property widget headers
    ${MOD_DIR}/qt/propertywidgets/colormappropertywidget.h
    ${MOD_DIR}/qt/propertywidgets/plotdatapropertywidget.h
    ${MOD_DIR}/qt/propertywidgets/plotentitiespropertywidget.h
    ${MOD_DIR}/qt/propertywidgets/plotpredicatepropertywidget.h
    ${MOD_DIR}/qt/propertywidgets/plotselectionpropertywidget.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/processorwidgets/plotprocessorwidgetfactory.h
    ${MOD_DIR}/qt/propertywidgets/plotpropertywidgetfactory.h
)
