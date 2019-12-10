/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <iostream>

#include "tgt/filesystem.h"
#include "tgt/logmanager.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/network/processornetwork.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE processornetworktest.cpp
#include <boost/test/unit_test.hpp>

using namespace voreen;
using std::vector;
using tgt::FileSystem;

//
// Check Functions
//
bool inVector(const vector<Processor*>& vec, const Processor* proc) {
    vector<Processor*>::const_iterator pos = std::find(vec.begin(), vec.end(), proc);
    return (pos != vec.end());
}

bool beforeInVector(const vector<Processor*>& vec, const Processor* first, const Processor* second) {
    vector<Processor*>::const_iterator firstPos = std::find(vec.begin(), vec.end(), first);
    vector<Processor*>::const_iterator secondPos = std::find(vec.begin(), vec.end(), second);
    return (secondPos < vec.end()) && (firstPos < secondPos);
}

bool isSizeLink(PropertyLink* link) {
    return (dynamic_cast<LinkEvaluatorRenderSize*>(link->getLinkEvaluator()) != 0);
}

int numSizeLinks(const std::vector<PropertyLink*> links) {
    int result = 0;
    for (size_t i=0; i<links.size(); i++)
        if (isSizeLink(links.at(i)))
            result++;
    return result;
}

bool portsSizeLinked(RenderPort* outport, RenderPort* inport) {
    tgtAssert(outport && inport, "null pointer passed");
    tgtAssert(outport->isOutport(), "passed port is not an outport");
    tgtAssert(inport->isInport(), "passed port is not an inport");

    RenderSizeReceiveProperty* receiveProp = outport->getSizeReceiveProperty();
    RenderSizeOriginProperty* originProp = inport->getSizeOriginProperty();
    if (!receiveProp || !originProp)
        return false;

    return (originProp->getLink(receiveProp) && isSizeLink(originProp->getLink(receiveProp)));
}

bool portsSizeLinked(RenderPort& outport, RenderPort& inport) {
    return portsSizeLinked(&outport, &inport);
}

class BasicPort : public Port {
public:
    BasicPort(PortDirection direction, const std::string& id, const std::string& guiName = "", bool allowMultipleConnections = false,
                Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT)
        : Port(direction, id, guiName, allowMultipleConnections, invalidationLevel)
    {}


    Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const{return new BasicPort(direction,id,guiName);}
    std::string getClassName() const { return "BasicPort";}
    bool hasData() const {return false;}
};

//
// Processor Mocks
//

/**
 * Default behaviour of MockProcessor:
 * - is invalid
 * - is ready
 * - invalidates its outports on process()
 */
class MockProcessor : public RenderProcessor {
public:
    MockProcessor()
        : RenderProcessor()
    {
        invalidate();
        ready = true;
    }

    virtual bool isReady() const {
        return ready;
    }

    virtual void process() {
        const std::vector<Port*> outports = getOutports();
        for (size_t i=0; i<outports.size(); i++)
            outports.at(i)->invalidatePort();
    }

    // workaround for linker error with precompiled headers and GCC 5
    virtual std::string getProgressMessage() const {
        if (!progressBars_.empty())
            return progressBars_.front()->getProgressMessage();
        else
            return "";
    }

    std::string getCategory() const  { return "Cat"; }
    void setDescriptions() {}

    bool ready;
};

/// Processor with active (size receiver) and passive render outport
class StartProcessor : public MockProcessor {
public:
    StartProcessor() :
        MockProcessor(),
        outportActive_(Port::OUTPORT, "outportActive", "outportActive", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
        outportPassive_(Port::OUTPORT, "outportPassive")
    {
        addPort(outportActive_);
        addPort(outportPassive_);
    }
    Processor* create() const        { return new StartProcessor(); }
    std::string getClassName() const { return "StartProcessor";     }

    RenderPort outportActive_;
    RenderPort outportPassive_;
};

/// Processor with two render inports and two render outports, one active/one passive each
class MiddleProcessor : public MockProcessor {
public:
    MiddleProcessor() :
        MockProcessor(),
        inportActive_(Port::INPORT, "inportActive", "inportActive", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
        inportActive2_(Port::INPORT, "inportActive2", "inportActive2", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
        inportPassive_(Port::INPORT, "inportPassive"),
        inportNonRender_(Port::INPORT, "inportNonRender"),
        outportActive_(Port::OUTPORT, "outportActive", "outportActive", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
        outportActive2_(Port::OUTPORT, "outportActive2", "outportActive2", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
        outportPassive_(Port::OUTPORT, "outportPassive"),
        outportNonRender_(Port::OUTPORT, "outportNonRender")
    {
        addPort(inportActive_);
        addPort(inportActive2_);
        addPort(inportPassive_);
        addPort(inportNonRender_);

        addPort(outportActive_);
        addPort(outportActive2_);
        addPort(outportPassive_);
        addPort(outportNonRender_);
    }
    Processor* create() const        { return new MiddleProcessor(); }
    std::string getClassName() const { return "MiddleProcessor";     }

    RenderPort inportActive_;
    RenderPort inportActive2_;
    RenderPort inportPassive_;

    RenderPort outportActive_;
    RenderPort outportActive2_;
    RenderPort outportPassive_;

    BasicPort inportNonRender_;
    BasicPort outportNonRender_;
};

/// Processor with two render inports (one active, one passive), marked as end processor
class EndProcessor : public MockProcessor {
public:
    EndProcessor() :
        MockProcessor(),
        inportActive_(Port::INPORT, "inportActive", "inportActive", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
        inportPassive_(Port::INPORT, "inportPassive")
    {
       addPort(inportActive_);
       addPort(inportPassive_);
    }
    Processor* create() const        { return new EndProcessor(); }
    std::string getClassName() const { return "EndProcessor";     }

    RenderPort inportActive_;
    RenderPort inportPassive_;
};

/// Processor with two normal render inports, two loop inports, two render outports, one active/passive each
class LoopInitiator : public MockProcessor {
public:
    LoopInitiator() :
        MockProcessor(),
        inportActive_(Port::INPORT, "inportActive", "inportActive", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
        inportPassive_(Port::INPORT, "inportPassive", "inportPassive"),
        loopInportActive_(Port::INPORT, "loopInportActive", "loopInportActive", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
        loopInportPassive_(Port::INPORT, "loopInportPassive", "loopInportPassive"),
        outportActive_(Port::OUTPORT, "outportActive", "outportActive", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
        outportPassive_(Port::OUTPORT, "outportPassive", "outportPassive")
    {
        loopInportActive_.setLoopPort(true);
        loopInportPassive_.setLoopPort(true);
        addPort(inportActive_);
        addPort(inportPassive_);
        addPort(loopInportActive_);
        addPort(loopInportPassive_);
        addPort(outportActive_);
        addPort(outportPassive_);
    }
    Processor* create() const        { return new LoopInitiator(); }
    std::string getClassName() const { return "LoopInitiator";     }

    RenderPort inportActive_;
    RenderPort inportPassive_;

    RenderPort loopInportActive_;
    RenderPort loopInportPassive_;

    RenderPort outportActive_;
    RenderPort outportPassive_;
};

/// Processor with two render inports, two normal render outports, two render loop outports, one active/passive each
class LoopFinalizer : public MockProcessor {
public:
    LoopFinalizer() :
      MockProcessor(),
          inportActive_(Port::INPORT, "inportActive", "inportActive", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
          inportPassive_(Port::INPORT, "inportPassive", "inportPassive"),
          outportActive_(Port::OUTPORT, "outportActive", "outportActive", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
          outportPassive_(Port::OUTPORT, "outportPassive", "outportPassive"),
          loopOutportActive_(Port::OUTPORT, "loopOutportActive", "loopOutportActive", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
          loopOutportPassive_(Port::OUTPORT, "loopOutportPassive", "loopOutportPassive")
      {
          loopOutportActive_.setLoopPort(true);
          loopOutportPassive_.setLoopPort(true);
          addPort(inportActive_);
          addPort(inportPassive_);
          addPort(outportActive_);
          addPort(outportPassive_);
          addPort(loopOutportActive_);
          addPort(loopOutportPassive_);
      }
      Processor* create() const        { return new LoopFinalizer(); }
      std::string getClassName() const { return "LoopFinalizer";     }

      RenderPort inportActive_;
      RenderPort inportPassive_;

      RenderPort outportActive_;
      RenderPort outportPassive_;

      RenderPort loopOutportActive_;
      RenderPort loopOutportPassive_;
};


//
// Test Fixtures
//

// Global setup & tear down
struct GlobalFixture {
    GlobalFixture() {
        tgt::Singleton<tgt::LogManager>::init();
        new VoreenApplication("networkevaluatortest", "networkevaluatortest", "NetworkEvaluatorTest",
            0, 0, VoreenApplication::APP_NONE);
        //VoreenApplication::app()->initialize();
    }
    ~GlobalFixture() {
        //VoreenApplication::app()->deinitialize();
        delete VoreenApplication::app();
        tgt::Singleton<tgt::LogManager>::deinit();
    }
};

/**
 * Executed before/after each test case:
 * - constructs two instances of each processor mock
 * - creates the ProcessorNetwork and adds the processor mocks
 * - creates the NetworkEvaluator and registers the EvaluationOrderRecorder, without assigning the network
 * - frees all created resources after test case execution
 */
struct TestFixture {
    StartProcessor *startProcessor, *startProcessor2;
    MiddleProcessor *middleProcessor, *middleProcessor2;
    EndProcessor *endProcessor, *endProcessor2;

    LoopInitiator *loopInitiator, *loopInitiator2;
    LoopFinalizer *loopFinalizer, *loopFinalizer2;

    ProcessorNetwork* network;

    TestFixture() {
        startProcessor = new StartProcessor();
        startProcessor2 = new StartProcessor();
        middleProcessor = new MiddleProcessor();
        middleProcessor2 = new MiddleProcessor();
        endProcessor = new EndProcessor();
        endProcessor2 = new EndProcessor();

        loopInitiator = new LoopInitiator();
        loopInitiator2 = new LoopInitiator();
        loopFinalizer = new LoopFinalizer();
        loopFinalizer2 = new LoopFinalizer();

        network = new ProcessorNetwork();

        network->addProcessor(startProcessor, "StartProcessor");
        network->addProcessor(startProcessor2, "StartProcessor2");
        network->addProcessor(middleProcessor, "MiddleProcessor");
        network->addProcessor(middleProcessor2, "MiddleProcessor2");
        network->addProcessor(endProcessor, "EndProcessor");
        network->addProcessor(endProcessor2, "EndProcessor2");

        network->addProcessor(loopInitiator, "LoopInitiator");
        network->addProcessor(loopInitiator2, "LoopInitiator2");
        network->addProcessor(loopFinalizer, "LoopFinalizer");
        network->addProcessor(loopFinalizer2, "LoopFinalizer2");
    }

    ~TestFixture() {
        delete network;
        // processor network deletes processors
    }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

//
// Test Cases
//
BOOST_FIXTURE_TEST_SUITE(SizeLinkingTests, TestFixture)

/****************************************************************************************/
/* --- Network notation ---                                                             */
/*                                                                                      */
/* - StartProcessor|a->a|MiddleProcessor:                                               */
/*      active (size receiving) render outport of StartProcessor is connected to        */
/*      active (size origin) render inport of MiddleProcessor                           */
/* - StartProcessor|p->p|MiddleProcessor:                                               */
/*      passive (non size receiving) render outport of StartProcessor is connected to   */
/*      passive (non size origin) render inport of MiddleProcessor                      */
/* - StartProcessor|a->p|MiddleProcessor, StartProcessor|p->a|MiddleProcessor:          */
/*      accordingly                                                                     */
/*                                                                                      */
/****************************************************************************************/

// Network = (StartProcessor)
BOOST_AUTO_TEST_CASE(singleProcessor)
{
    vector<Processor*> processors;
    processors.push_back(startProcessor);

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(processors) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor,MiddleProcessor,EndProcessor)
BOOST_AUTO_TEST_CASE(unconnectedGraph)
{
    vector<Processor*> processors;
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor);
    processors.push_back(endProcessor);

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(processors) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a,a|MiddleProcessor|a,a|EndProcessor)
//                           ^-/                 ^-/
BOOST_AUTO_TEST_CASE(unconnectedGraph_RemoveLinksFromNetwork)
{
    // create size links
    network->createPropertyLink(middleProcessor->inportActive_.getSizeOriginProperty(),
        startProcessor->outportActive_.getSizeReceiveProperty(),
        new LinkEvaluatorRenderSize());
    network->createPropertyLink(endProcessor->inportActive_.getSizeOriginProperty(),
        middleProcessor->outportActive_.getSizeReceiveProperty(),
        new LinkEvaluatorRenderSize());

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromSubNetwork(network->getProcessors()) == 2);

    // check that links have been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a,a|MiddleProcessor|a,a|EndProcessor)
//                           ^-/                 ^-/
BOOST_AUTO_TEST_CASE(unconnectedGraph_RemoveLinksOverConnection)
{
    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&startProcessor->outportActive_, &middleProcessor->inportActive_) == 1);

    // check that link spanning the startProcssor->middleProcessor connection has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a,a|MiddleProcessor|a,a|EndProcessor)
//                           ^-/                 ^-/
BOOST_AUTO_TEST_CASE(unconnectedGraph_RemoveLinksOverFalseConnection)
{
    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&startProcessor->outportActive_, &endProcessor->inportActive_) == 0);

    // no size link spans the (imaginary) connection from startProcessor to endProcessor
    // => no size link should have been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(twoProcessorsActive)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|EndProcessor)
//                           ^--/
BOOST_AUTO_TEST_CASE(twoProcessorsActive_ExistingLink)
{
    // create port connections and size links
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors(), false) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|EndProcessor)
//                           ^--/
BOOST_AUTO_TEST_CASE(twoProcessorsActive_ExistingLinkReplace)
{
    // create port connections and size links
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors(), true) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(twoProcessorsActive_OverConnection)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&startProcessor->outportActive_, &endProcessor->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|EndProcessor)
//                           ^--/
BOOST_AUTO_TEST_CASE(twoProcessorsActive_RemoveLinksFromNetwork)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromSubNetwork(network->getProcessors()) == 1);

    // check that links have been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->a|EndProcessor)
//                           ^--/
BOOST_AUTO_TEST_CASE(twoProcessorsActive_RemoveLinksOverConnection)
{
    vector<Processor*> processors;
    processors.push_back(startProcessor);
    processors.push_back(endProcessor);

    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&startProcessor->outportActive_, &endProcessor->inportActive_) == 1);

    // check that links have been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|p->p|EndProcessor)
BOOST_AUTO_TEST_CASE(twoProcessorsPassive)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportPassive_, endProcessor->inportPassive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|p->p|EndProcessor)
BOOST_AUTO_TEST_CASE(twoProcessorsPassive_RemoveLinksFromSubNetwork)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportPassive_, endProcessor->inportPassive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromSubNetwork(network->getProcessors()) == 0);
}

// Network = (StartProcessor|p->p|EndProcessor)
BOOST_AUTO_TEST_CASE(twoProcessorsPassive_RemoveLinksOverConnection)
{
    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&startProcessor->outportPassive_, &endProcessor->inportPassive_) == 0);
}

// Network = (StartProcessor|a->p|EndProcessor)
BOOST_AUTO_TEST_CASE(twoProcessorsActiveToPassive)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportPassive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|p->a|EndProcessor)
BOOST_AUTO_TEST_CASE(twoProcessorsPassiveToActive)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportPassive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->a|MiddleProcessor|a->a|MiddleProcessor2|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineFullActive)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 3);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 3);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a->a|MiddleProcessor2|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_OverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&middleProcessor->outportActive_, &middleProcessor2->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a->a|MiddleProcessor2|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_CreateLinksForProcessor)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create render size link
    BOOST_CHECK(network->createRenderSizeLinksForProcessor(middleProcessor) == 2);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
}

// Network = (StartProcessor|a->a|~MiddleProcessor~|a->a|MiddleProcessor2|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_SkippedMiddle)
{
    vector<Processor*> processors;
    processors.push_back(startProcessor);
    //processors.push_back(middleProcessor);  //<<!!!
    processors.push_back(middleProcessor2);
    processors.push_back(endProcessor);

    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(processors) == 2);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|~MiddleProcessor~|a->a|~MiddleProcessor2~|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_SkippedMiddle2)
{
    vector<Processor*> processors;
    processors.push_back(startProcessor);
    //processors.push_back(middleProcessor);  //<<!!!
    //processors.push_back(middleProcessor2); //<<!!!
    processors.push_back(endProcessor);

    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(processors) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a->a|MiddleProcessor2|a->a|EndProcessor)
//                                                ^------------------------/
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_ReplaceExisting)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size link
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor->outportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors(), true) == 3);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 3);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a->a|MiddleProcessor2|a->a|EndProcessor)
//                                                ^------------------------/
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_NotReplaceExisting)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size link
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor->outportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors(), false) == 2);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 3);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a->a|MiddleProcessor2|a->a|EndProcessor)
//                           ^--/                 ^--/                  ^--/
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_RemoveLinksFromNetwork)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor2->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromSubNetwork(network->getProcessors()) == 3);

    // check that links have been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (~StartProcessor~|a->a|MiddleProcessor|a->a|MiddleProcessor2|a->a|~EndProcessor~)
//                             ^--/                 ^--/                  ^--/
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_RemoveMiddleLinksFromIncompleteNetwork)
{
    vector<Processor*> processors;
    //processors.push_back(startProcessor);  //<<!!!
    processors.push_back(middleProcessor);
    processors.push_back(middleProcessor2);
    //processors.push_back(endProcessor);    //<<!!!

    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor2->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromSubNetwork(processors) == 1);

    // check that link has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (~StartProcessor~|a->a|MiddleProcessor|a->a|~MiddleProcessor2~|a->a|EndProcessor)
//                             ^--/                 ^--/                    ^--/
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_RemoveMiddleLinksFromIncompleteNetwork2)
{
    vector<Processor*> processors;
    //processors.push_back(startProcessor);   //<<!!!
    processors.push_back(middleProcessor);
    //processors.push_back(middleProcessor2); //<<!!!
    processors.push_back(endProcessor);

    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor2->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromSubNetwork(processors) == 0);

    // check that links are still present
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 3);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a->a|MiddleProcessor2|a->a|EndProcessor)
//                           ^--/                 ^--/                  ^--/
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_RemoveLinksOverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor2->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&middleProcessor->outportActive_, &middleProcessor2->inportActive_) == 1);

    // check that (only) link spanning the middleProcessor->middleProcessor2 connection has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a->a|MiddleProcessor2|a->a|EndProcessor)
//                           ^--/                 ^--/                  ^--/
BOOST_AUTO_TEST_CASE(linearPipelineFullActive_RemoveLinksFromProcessor)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor2->outportActive_));

    // remove size link from processor
    BOOST_CHECK(network->removeRenderSizeLinksFromProcessor(middleProcessor) == 2);

    // check that links have been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|p->p|MiddleProcessor|p->p|MiddleProcessor2|p->p|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineFullPassive)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportPassive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportPassive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|p->p|MiddleProcessor|p->p|MiddleProcessor2|p->p|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineFullPassive_CreateLinksForProcessor)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportPassive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportPassive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksForProcessor(middleProcessor) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|p->p|MiddleProcessor|p->p|MiddleProcessor2|p->p|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineFullPassive_RemoveLinksFromProcessor)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportPassive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportPassive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromProcessor(middleProcessor) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive_OverPassiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&middleProcessor->outportPassive_, &middleProcessor2->inportPassive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive_OverActiveToPassiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&startProcessor->outportActive_, &middleProcessor->inportPassive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive_OverPassiveToActiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&middleProcessor2->outportPassive_, &endProcessor->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
//                           ^---------------------------------------------/
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive_ReplaceExisting)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // configure existing size link
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors(), true) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
//                           ^---------------------------------------------/
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive_NotReplaceExisting)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // configure existing size link
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors(), false) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
//                           ^---------------------------------------------/
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive_RemoveLinksFromSubNetwork)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromSubNetwork(network->getProcessors()) == 1);

    // check that link has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
//                           ^---------------------------------------------/
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive_RemoveLinksOverPassiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // configure existing size link
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&middleProcessor->outportPassive_, &middleProcessor2->inportPassive_) == 1);

    // check that link has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
//                           ^---------------------------------------------/
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive_RemoveLinksOverActiveToPassiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&startProcessor->outportActive_, &middleProcessor->inportPassive_) == 1);

    // check that link has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->p|MiddleProcessor|p->p|MiddleProcessor2|p->a|EndProcessor)
//                           ^---------------------------------------------/
BOOST_AUTO_TEST_CASE(linearPipelineMiddlePassive_RemoveLinksOverPassiveToActiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&middleProcessor2->outportPassive_, &endProcessor->inportActive_) == 1);

    // check that link has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|p->p|MiddleProcessor|a->a|MiddleProcessor2|p->p|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineMiddleActive)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportPassive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
}

// Network = (StartProcessor|p->p|MiddleProcessor|a->a|MiddleProcessor2|p->p|EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipelineMiddleActive_OverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportPassive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportPassive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&middleProcessor->outportActive_, &middleProcessor2->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
}

// Network = (StartProcessor|a->a|EndProcessor )
//                           \->a|EndProcessor2)
BOOST_AUTO_TEST_CASE(multipleOrigins)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor2->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 0);

    // check results: no size links created due to ambiguity
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->a|EndProcessor )
//                           \->a|EndProcessor2)
BOOST_AUTO_TEST_CASE(multipleOrigins_OverConnection)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor2->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&startProcessor->outportActive_, &endProcessor->inportActive_) == 0);

    // check results: no size links created due to ambiguity
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->a|EndProcessor )
//                           \->a|EndProcessor2)
BOOST_AUTO_TEST_CASE(multipleOrigins_RemoveLinksOverConnection)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor2->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_,  startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor2->inportActive_, startProcessor->outportActive_));

    // remove size link
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&startProcessor->outportActive_, &endProcessor->inportActive_) == 1);

    // check that (only) size link startProcessor->endProcessor has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor2->inportActive_));
}

// Network = (StartProcessor|a-->a|EndProcessor)
//            StartProcessor2|a--^             )
BOOST_AUTO_TEST_CASE(multipleReceivers)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 2);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(startProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a-->a|EndProcessor)
//            StartProcessor2|a--^             )
//                            ^--/             )
BOOST_AUTO_TEST_CASE(multipleReceivers_ReplaceOneExisting)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size link
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor2->outportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors(), true) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(startProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a-->a|EndProcessor)
//            StartProcessor2|a--^             )
//                            ^--/             )
BOOST_AUTO_TEST_CASE(multipleReceivers_NotReplaceOneExisting)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size link
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor2->outportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors(), false) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(startProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a-->a|EndProcessor)
//            StartProcessor2|a--^             )
BOOST_AUTO_TEST_CASE(multipleReceivers_OverConnection)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&startProcessor->outportActive_, &endProcessor->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a-->a|EndProcessor)
//            StartProcessor2|a--^             )
BOOST_AUTO_TEST_CASE(multipleReceivers_CreateLinksForEndProcessor)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksForProcessor(endProcessor) == 2);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(startProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a-->a|EndProcessor)
//            StartProcessor2|a--^             )
BOOST_AUTO_TEST_CASE(multipleReceivers_CreateLinksForStartProcessor)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksForProcessor(startProcessor) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a-->a|EndProcessor)
//            StartProcessor2|a--^             )
BOOST_AUTO_TEST_CASE(multipleReceivers_RemoveLinksOverConnection)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor2->outportActive_));

    // remove size link
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&startProcessor->outportActive_, &endProcessor->inportActive_) == 1);

    // check that (only) size link startProcessor->endProcessor has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a-->a|EndProcessor)
//            StartProcessor2|a--^             )
BOOST_AUTO_TEST_CASE(multipleReceivers_RemoveLinksFromStartProcessor)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor2->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromProcessor(startProcessor) == 1);

    // check that (only) size link startProcessor->endProcessor has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a-->a|EndProcessor)
//            StartProcessor2|a--^             )
BOOST_AUTO_TEST_CASE(multipleReceivers_RemoveLinksFromEndProcessor)
{
    // create port connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor2->outportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromProcessor(endProcessor) == 2);

    // check that both links have been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->p|MiddleProcessor|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(receiverBetweenReceiverAndOrigin)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|MiddleProcessor|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(receiverBetweenReceiverAndOrigin_OverActiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&middleProcessor->outportActive_, &endProcessor->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|MiddleProcessor|a->a|EndProcessor)
BOOST_AUTO_TEST_CASE(receiverBetweenReceiverAndOrigin_OverActiveToPassiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&startProcessor->outportActive_, &middleProcessor->inportPassive_) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->a|MiddleProcessor|a1|a2->a2|a1|MiddleProcessor2|a->a|EndProcessor)
//                                                 \--------^
BOOST_AUTO_TEST_CASE(dualConnectionActive)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive2_, middleProcessor2->inportActive2_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 4);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 4);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive2_, middleProcessor2->inportActive2_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a1|a2->a2|a1|MiddleProcessor2|a->a|EndProcessor)
//                                                 \--------^
BOOST_AUTO_TEST_CASE(dualConnectionActive_OverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive2_, middleProcessor2->inportActive2_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&middleProcessor->outportActive_, &middleProcessor2->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a1|a2->a2|a1|MiddleProcessor2|a->a|EndProcessor)
//                                                 \--------^
BOOST_AUTO_TEST_CASE(dualConnectionActive_CreateLinksForProcessor)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive2_, middleProcessor2->inportActive2_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksForProcessor(middleProcessor) == 3);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 3);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive2_, middleProcessor2->inportActive2_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a1|a2->a2|a1|MiddleProcessor2|a->a|EndProcessor)
//                                                 \--------^
BOOST_AUTO_TEST_CASE(dualConnectionActive_RemoveLinksOverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive2_, middleProcessor2->inportActive2_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive2_, middleProcessor->outportActive2_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor2->outportActive_));

    // remove link
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&middleProcessor->outportActive_, &middleProcessor2->inportActive_) == 1);

    // check that (only) link spanning the connection has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 3);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    //BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive2_, middleProcessor2->inportActive2_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a1|a2->a2|a1|MiddleProcessor2|a->a|EndProcessor)
//                                                 \--------^
BOOST_AUTO_TEST_CASE(dualConnectionActive_RemoveLinksFromProcessor)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive2_, middleProcessor2->inportActive2_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive2_, middleProcessor->outportActive2_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor2->outportActive_));

    // remove links
    BOOST_CHECK(network->removeRenderSizeLinksFromProcessor(middleProcessor) == 3);

    // check that (only) links from/to middle processor have been deleted
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a|p->p|a|MiddleProcessor2|a->a|EndProcessor)
//                                                 \-----^
BOOST_AUTO_TEST_CASE(dualConnectionOnePassive)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 3);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 3);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    //BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive2_, middleProcessor2->inportActive2_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a|p->p|a|MiddleProcessor2|a->a|EndProcessor)
//                                                 \-----^
BOOST_AUTO_TEST_CASE(dualConnectionOnePassive_OverActiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&middleProcessor->outportActive_, &middleProcessor2->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a|p->p|a|MiddleProcessor2|a->a|EndProcessor)
//                                                 \-----^
BOOST_AUTO_TEST_CASE(dualConnectionOnePassive_OverPassiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&middleProcessor->outportPassive_, &middleProcessor2->inportPassive_) == 0);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->a|MiddleProcessor|a|p->p|a|MiddleProcessor2|a->a|EndProcessor)
//                                                 \-----^
BOOST_AUTO_TEST_CASE(dualConnectionOnePassive_RemoveLinksOverActiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor2->outportActive_));

    // remove link
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&middleProcessor->outportActive_, &middleProcessor2->inportActive_) == 1);

    // check that (only) link spanning the connection has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|MiddleProcessor|a|p->p|a|MiddleProcessor2|a->a|EndProcessor)
//                                                 \-----^
BOOST_AUTO_TEST_CASE(dualConnectionOnePassive_RemoveLinksOverPassiveConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, middleProcessor2->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outportActive_, endProcessor->inportActive_));

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor2->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, middleProcessor2->outportActive_));

    // remove link
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&middleProcessor->outportPassive_, &middleProcessor2->inportPassive_) == 0);

    // check no link has been removed
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 3);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, middleProcessor2->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor2->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|(la|LoopInitiator|a->a|MiddleProcessor|a->a|LoopFinalizer|la)|a->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(activeLoop_CreateSizeLinksInSubnetwork)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportActive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, loopFinalizer->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportActive_, loopInitiator->loopInportActive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 4);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 4);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, loopInitiator->inportActive_));
    BOOST_CHECK(portsSizeLinked(loopInitiator->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, loopFinalizer->inportActive_));
    BOOST_CHECK(portsSizeLinked(loopFinalizer->outportActive_, endProcessor->inportActive_));

    // NOTE: loop port connections are expected to be ignored, even if they are active!
    BOOST_CHECK(!portsSizeLinked(loopFinalizer->loopOutportActive_, loopInitiator->loopInportActive_));
}

// Network = (StartProcessor|a->a|(la|LoopInitiator|a->a|MiddleProcessor|a->a|LoopFinalizer|la)|a->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(activeLoop_RemoveSizeLinksInSubnetwork)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportActive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, loopFinalizer->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportActive_, loopInitiator->loopInportActive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(loopInitiator->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, loopInitiator->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(loopFinalizer->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, loopFinalizer->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(loopInitiator->loopInportActive_, loopFinalizer->loopOutportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksFromSubNetwork(network->getProcessors()) == 5);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->a|(la|LoopInitiator|a->a|MiddleProcessor|a->a|LoopFinalizer|la)|a->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(activeLoop_CreateSizeLinksOverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportActive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, loopFinalizer->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportActive_, loopInitiator->loopInportActive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&loopInitiator->outportActive_, &middleProcessor->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(loopInitiator->outportActive_, middleProcessor->inportActive_));
}

// Network = (StartProcessor|a->a|(la|LoopInitiator|a->a|MiddleProcessor|a->a|LoopFinalizer|la)|a->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(activeLoop_RemoveSizeLinksOverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportActive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, loopFinalizer->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportActive_, loopInitiator->loopInportActive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(loopInitiator->inportActive_, startProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(middleProcessor->inportActive_, loopInitiator->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(loopFinalizer->inportActive_, middleProcessor->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, loopFinalizer->outportActive_));
    BOOST_CHECK(network->createRenderSizeLink(loopInitiator->loopInportActive_, loopFinalizer->loopOutportActive_));

    // remove size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&loopInitiator->outportActive_, &middleProcessor->inportActive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 4);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, loopInitiator->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, loopFinalizer->inportActive_));
    BOOST_CHECK(portsSizeLinked(loopFinalizer->outportActive_, endProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(loopFinalizer->loopOutportActive_, loopInitiator->loopInportActive_));
}


// Network = (StartProcessor|a->p|(lp|LoopInitiator|p->p|MiddleProcessor|p->p|LoopFinalizer|lp)|p->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(passiveLoop_CreateSizeLinksInSubnetwork)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportPassive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, loopFinalizer->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportPassive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportPassive_, loopInitiator->loopInportPassive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|(lp|LoopInitiator|p->p|MiddleProcessor|p->p|LoopFinalizer|lp)|p->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(passiveLoop_RemoveSizeLinksFromSubnetwork)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportPassive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, loopFinalizer->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportPassive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportPassive_, loopInitiator->loopInportPassive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // create size links
    BOOST_CHECK(network->removeRenderSizeLinksFromSubNetwork(network->getProcessors()) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->p|(lp|LoopInitiator|p->p|MiddleProcessor|p->p|LoopFinalizer|lp)|p->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(passiveLoop_CreateSizeLinksOverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportPassive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, loopFinalizer->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportPassive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportPassive_, loopInitiator->loopInportPassive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&loopInitiator->outportPassive_, &middleProcessor->inportPassive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|(lp|LoopInitiator|p->p|MiddleProcessor|p->p|LoopFinalizer|lp)|p->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(passiveLoop_RemoveSizeLinksOverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportPassive_, middleProcessor->inportPassive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportPassive_, loopFinalizer->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportPassive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportPassive_, loopInitiator->loopInportPassive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // configure existing size links
    BOOST_CHECK(network->createRenderSizeLink(endProcessor->inportActive_, startProcessor->outportActive_));

    // create size links
    BOOST_CHECK(network->removeRenderSizeLinksOverConnection(&loopInitiator->outportPassive_, &middleProcessor->inportPassive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 0);
}

// Network = (StartProcessor|a->p|(lp|LoopInitiator|p->a|MiddleProcessor|a->p|LoopFinalizer|lp)|p->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(middleActiveLoop_CreateSizeLinksInSubnetwork)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportPassive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, loopFinalizer->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportPassive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportPassive_, loopInitiator->loopInportPassive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksWithinSubNetwork(network->getProcessors()) == 2);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 2);
    BOOST_CHECK(portsSizeLinked(startProcessor->outportActive_, middleProcessor->inportActive_));
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, endProcessor->inportActive_));
}

// Network = (StartProcessor|a->p|(lp|LoopInitiator|p->a|MiddleProcessor|a->p|LoopFinalizer|lp)|p->a|EndProcessor)
//                                  ^-------------------------x3----------------------------/
BOOST_AUTO_TEST_CASE(middleActiveLoop_CreateSizeLinksOverConnection)
{
    // configure network connections
    BOOST_CHECK(network->connectPorts(startProcessor->outportActive_, loopInitiator->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outportPassive_, middleProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outportActive_, loopFinalizer->inportPassive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outportPassive_, endProcessor->inportActive_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutportPassive_, loopInitiator->loopInportPassive_));
    loopInitiator->loopInportActive_.setNumLoopIterations(3);

    // create size links
    BOOST_CHECK(network->createRenderSizeLinksOverConnection(&middleProcessor->outportActive_, &loopFinalizer->inportPassive_) == 1);

    // check results
    BOOST_CHECK(numSizeLinks(network->getPropertyLinks()) == 1);
    BOOST_CHECK(portsSizeLinked(middleProcessor->outportActive_, endProcessor->inportActive_));
}
BOOST_AUTO_TEST_SUITE_END()
