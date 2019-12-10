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
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/port.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/networkgraph.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE networkevaluatortest.cpp
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
class MockProcessor : public Processor {
public:
    MockProcessor()
        : Processor()
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

/// Processor with two outports
class StartProcessor : public MockProcessor {
public:
    StartProcessor() :
        MockProcessor(),
        outport_(Port::OUTPORT, "outport"),
        outport2_(Port::OUTPORT, "outport2")
    {
        addPort(outport_);
        addPort(outport2_);
    }
    Processor* create() const        { return new StartProcessor(); }
    std::string getClassName() const { return "StartProcessor";     }

    BasicPort outport_;
    BasicPort outport2_;
};

/// Processor with two inports and two outports
class MiddleProcessor : public MockProcessor {
public:
    MiddleProcessor() :
        MockProcessor(),
        inport_(Port::INPORT, "inport"),
        inport2_(Port::INPORT, "inport2"),
        outport_(Port::OUTPORT, "outport"),
        outport2_(Port::OUTPORT, "outport2")
    {
        addPort(inport_);
        addPort(inport2_);
        addPort(outport_);
        addPort(outport2_);
    }
    Processor* create() const        { return new MiddleProcessor(); }
    std::string getClassName() const { return "MiddleProcessor";     }

    BasicPort inport_;
    BasicPort inport2_;

    BasicPort outport_;
    BasicPort outport2_;
};

/// Processor with two inports, marked as end processor
class EndProcessor : public MockProcessor {
public:
    EndProcessor() :
        MockProcessor(),
        inport_(Port::INPORT, "inport"),
        inport2_(Port::INPORT, "inport2")
    {
       addPort(inport_);
       addPort(inport2_);
    }
    Processor* create() const        { return new EndProcessor(); }
    std::string getClassName() const { return "EndProcessor";     }

    BasicPort inport_;
    BasicPort inport2_;
};

/// Processor with one normal inport, one loop inport, and one outport
class LoopInitiator : public MockProcessor {
public:
    LoopInitiator() :
        MockProcessor(),
        inport_(Port::INPORT, "inport"),
        loopInport_(Port::INPORT, "loopInport"),
        outport_(Port::OUTPORT, "outport")
    {
        loopInport_.setLoopPort(true);
        addPort(inport_);
        addPort(loopInport_);
        addPort(outport_);
    }
    Processor* create() const        { return new LoopInitiator(); }
    std::string getClassName() const { return "LoopInitiator";     }

    BasicPort inport_;
    BasicPort loopInport_;
    BasicPort outport_;
};

/// Processor with one normal inport, one normal outport, and one loop outport
class LoopFinalizer : public MockProcessor {
public:
    LoopFinalizer() :
        MockProcessor(),
        inport_(Port::INPORT, "inport"),
        outport_(Port::OUTPORT, "outport"),
        loopOutport_(Port::OUTPORT, "loopOutport")
    {
        addPort(inport_);
        addPort(outport_);
        loopOutport_.setLoopPort(true);
        addPort(loopOutport_);
    }
    Processor* create() const        { return new LoopFinalizer(); }
    std::string getClassName() const { return "LoopFinalizer";     }

    BasicPort inport_;
    BasicPort outport_;
    BasicPort loopOutport_;
};


//
// Test Fixtures
//

/// Is registered at the NetworkEvaluator and records the evaluation order.
class EvaluationOrderRecorder : public NetworkEvaluatorObserver {
public:
    virtual void beforeProcess(Processor* p) {
        evalOrder.push_back(p);
        if (p->getPort("loopInport"))
            loopIterations.push_back(p->getPort("loopInport")->getLoopIteration());
        else if (p->getPort("loopOutport"))
            loopIterations.push_back(p->getPort("loopOutport")->getLoopIteration());
        else
            loopIterations.push_back(-1);
    }
    virtual void afterProcess(Processor* /*p*/) {}

    std::vector<Processor*> evalOrder;
    std::vector<int> loopIterations;
};

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
    LoopInitiator *loopInitiator, *loopInitiator2;
    LoopFinalizer *loopFinalizer, *loopFinalizer2;
    EndProcessor *endProcessor, *endProcessor2;

    ProcessorNetwork* network;
    NetworkEvaluator* evaluator;
    EvaluationOrderRecorder* evalOrderRec;

    TestFixture() {
        startProcessor = new StartProcessor();
        startProcessor2 = new StartProcessor();
        middleProcessor = new MiddleProcessor();
        middleProcessor2 = new MiddleProcessor();
        loopInitiator = new LoopInitiator();
        loopInitiator2 = new LoopInitiator();
        loopFinalizer = new LoopFinalizer();
        loopFinalizer2 = new LoopFinalizer();
        endProcessor = new EndProcessor();
        endProcessor2 = new EndProcessor();

        network = new ProcessorNetwork();
        evaluator = new NetworkEvaluator(false);
        evalOrderRec = new EvaluationOrderRecorder();
        evaluator->addObserver(evalOrderRec);

        network->addProcessor(startProcessor, "StartProcessor");
        network->addProcessor(startProcessor2, "StartProcessor2");
        network->addProcessor(middleProcessor, "MiddleProcessor");
        network->addProcessor(middleProcessor2, "MiddleProcessor2");
        network->addProcessor(loopInitiator, "LoopInitiator");
        network->addProcessor(loopInitiator2, "LoopInitiator2");
        network->addProcessor(loopFinalizer, "LoopFinalizer");
        network->addProcessor(loopFinalizer2, "LoopFinalizer2");
        network->addProcessor(endProcessor, "EndProcessor");
        network->addProcessor(endProcessor2, "EndProcessor2");
    }

    ~TestFixture() {
        evaluator->addObserver(evalOrderRec);
        evaluator->setProcessorNetwork(nullptr, true);
        delete evaluator;
        delete evalOrderRec;
        delete network;
        // processor network deletes processors
    }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

//
// Test Cases
//
BOOST_FIXTURE_TEST_SUITE(NetworkGraphTests, TestFixture)

BOOST_AUTO_TEST_CASE(emptyGraph)
{
    vector<Processor*> processors;
    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().empty());
    BOOST_CHECK(!graph.containsProcessor(0));
    BOOST_CHECK(!graph.containsProcessor(startProcessor));

    BOOST_CHECK(graph.getPredecessors(endProcessor).empty());
    BOOST_CHECK(graph.getSuccessors(startProcessor).empty());

    NetworkGraph* subGraph = graph.getSubGraph();
    BOOST_CHECK(subGraph->getProcessors().empty());
    delete subGraph;

    NetworkGraph* transposed = graph.getTransposed();
    BOOST_CHECK(transposed->getProcessors().empty());
    delete transposed;
}

// Network = (StartProcessor)
BOOST_AUTO_TEST_CASE(singleProcessor)
{
    vector<Processor*> processors;
    processors.push_back(startProcessor);
    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 1);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(!graph.containsProcessor(endProcessor));

    // processor predecessors/successors
    BOOST_CHECK(graph.getPredecessors(endProcessor).empty()); //< end processor not contained by network
    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor, startProcessor));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    // port successors
    BOOST_CHECK(graph.getSuccessors(&startProcessor->outport_).empty());
    BOOST_CHECK(graph.getSuccessors(&startProcessor->outport2_).empty());

    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport_, &startProcessor->outport2_));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport_, startProcessor));
    BOOST_CHECK(graph.isPathElement(&startProcessor->outport_, &startProcessor->outport_, &startProcessor->outport_));
    BOOST_CHECK(!graph.isPathElement(&startProcessor->outport_, &startProcessor->outport_, &startProcessor->outport2_));

}

// Network = (StartProcessor,MiddleProcessor,EndProcessor)
BOOST_AUTO_TEST_CASE(unconnectedGraph)
{
    vector<Processor*> processors;
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor);
    processors.push_back(endProcessor);
    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 3);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // unconnected graph => processors have no predecessors/successors
    BOOST_CHECK(graph.getPredecessors(endProcessor).size() == 1 && graph.getPredecessors(endProcessor).count(endProcessor));
    BOOST_CHECK(graph.getSuccessors(startProcessor).size() == 1 && graph.getSuccessors(startProcessor).count(startProcessor));
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport_).empty());
    BOOST_CHECK(graph.getSuccessors(&startProcessor->outport_).empty());
    BOOST_CHECK(!graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor, endProcessor));

    // no port successors/predecessors
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport_));

    // port successors/predecessors
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, &endProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport_));

    BOOST_CHECK(graph.isSuccessor(&middleProcessor->inport_, &middleProcessor->outport_));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, &middleProcessor->inport_));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->inport_, middleProcessor));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, middleProcessor));
    BOOST_CHECK(!graph.isSuccessor(middleProcessor, &middleProcessor->inport_));
    BOOST_CHECK(graph.isSuccessor(middleProcessor, &middleProcessor->outport_));
    BOOST_CHECK(graph.isPathElement(middleProcessor, &middleProcessor->inport_, &middleProcessor->outport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &middleProcessor->outport_, &middleProcessor->inport_));

}

// Network = (StartProcessor->MiddleProcessor->MiddleProcessor2->EndProcessor)
BOOST_AUTO_TEST_CASE(linearPipeline)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(middleProcessor2);
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, endProcessor->inport_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor2));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // processors are connected => processors are predecessors/successors of each other
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(middleProcessor)  &&
                graph.getPredecessors(endProcessor).count(middleProcessor2) &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor).count(middleProcessor2) &&
                graph.getSuccessors(startProcessor).count(endProcessor)
        );
    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));

    // port predecessors/successors
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport_).count(startProcessor)     &&
                graph.getPredecessors(&endProcessor->inport_).count(middleProcessor)    &&
                graph.getPredecessors(&endProcessor->inport_).count(middleProcessor2)   &&
                !graph.getPredecessors(&endProcessor->inport_).count(endProcessor)
        );
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport2_).empty());

    BOOST_CHECK(!graph.getSuccessors(&startProcessor->outport_).count(startProcessor)   &&
                graph.getSuccessors(&startProcessor->outport_).count(middleProcessor)   &&
                graph.getSuccessors(&startProcessor->outport_).count(middleProcessor2)  &&
                graph.getSuccessors(&startProcessor->outport_).count(endProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(&startProcessor->outport2_).empty());

    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, endProcessor));

    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, &middleProcessor->inport_));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, &middleProcessor->outport_));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, &middleProcessor->outport_));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, &middleProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, &middleProcessor->inport2_));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, &middleProcessor->inport_));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->inport_, endProcessor));
    BOOST_CHECK(graph.isSuccessor(middleProcessor, &middleProcessor2->outport_));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->inport_, &middleProcessor2->inport_));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->inport_, &middleProcessor2->inport2_));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->inport_, &middleProcessor2->outport_));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->inport_, &middleProcessor2->outport2_));

    BOOST_CHECK(graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport2_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport2_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport2_, &endProcessor->inport2_));
    BOOST_CHECK(!graph.isPathElement(startProcessor, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(endProcessor, &startProcessor->outport_, &endProcessor->inport_));

    BOOST_CHECK(graph.isPathElement(middleProcessor2, &middleProcessor->inport_, &middleProcessor2->outport_));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, &middleProcessor->inport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, &middleProcessor->outport2_, &middleProcessor2->outport_));

    BOOST_CHECK(graph.isPathElement(&middleProcessor->inport_, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(graph.isPathElement(&middleProcessor->outport_, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(&middleProcessor->inport2_, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(&middleProcessor->outport2_, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(&middleProcessor->inport_, &startProcessor->outport_, &endProcessor->inport2_));
    BOOST_CHECK(!graph.isPathElement(&middleProcessor->inport_, &startProcessor->outport2_, &endProcessor->inport_));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 4);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, middleProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor2, endProcessor));
}

// Network = (StartProcessor->EndProcessor,MiddleProcessor)
BOOST_AUTO_TEST_CASE(skip)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, endProcessor->inport_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 3);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // start/end processor are predecessors/successors or each other, but middleProcessor is unconnected
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(endProcessor)     &&
                !graph.getPredecessors(endProcessor).count(middleProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(endProcessor)     &&
                !graph.getSuccessors(startProcessor).count(middleProcessor)
        );
    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor, endProcessor));

    // port predecessors/successors
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport_).count(startProcessor)     &&
                !graph.getPredecessors(&endProcessor->inport_).count(endProcessor)      &&
                !graph.getPredecessors(&endProcessor->inport_).count(middleProcessor)
    );
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport2_).empty());
    BOOST_CHECK(graph.getPredecessors(&middleProcessor->inport_).empty());

    BOOST_CHECK(!graph.getSuccessors(&startProcessor->outport_).count(startProcessor)   &&
                graph.getSuccessors(&startProcessor->outport_).count(endProcessor)      &&
                !graph.getSuccessors(&startProcessor->outport_).count(middleProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(&startProcessor->outport2_).empty());
    BOOST_CHECK(graph.getSuccessors(&middleProcessor->outport_).empty());
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, endProcessor));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport_));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 3);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, endProcessor));
}

// Network = (StartProcessor->MiddleProcessor->EndProcessor )
//           (             \------------------>EndProcessor2)
BOOST_AUTO_TEST_CASE(halfSkip)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor2);
    processors.push_back(endProcessor);
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(startProcessor->outport2_, endProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, endProcessor->inport_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor2));

    // check predecessors/successors
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)    &&
                graph.getPredecessors(endProcessor).count(middleProcessor)   &&
                graph.getPredecessors(endProcessor).count(endProcessor)      &&
                !graph.getPredecessors(endProcessor).count(endProcessor2)
        );

    BOOST_CHECK(graph.getPredecessors(endProcessor2).count(startProcessor)   &&
                !graph.getPredecessors(endProcessor2).count(middleProcessor) &&
                graph.getPredecessors(endProcessor2).count(endProcessor2)    &&
                !graph.getPredecessors(endProcessor2).count(endProcessor)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)    &&
                graph.getSuccessors(startProcessor).count(middleProcessor)   &&
                graph.getSuccessors(startProcessor).count(endProcessor)      &&
                graph.getSuccessors(startProcessor).count(endProcessor2)
        );
    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor2));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor, endProcessor2));

    // port predecessors/successors
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport_).count(startProcessor)    &&
                graph.getPredecessors(&endProcessor->inport_).count(middleProcessor)   &&
                !graph.getPredecessors(&endProcessor->inport_).count(endProcessor)     &&
                !graph.getPredecessors(&endProcessor->inport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport2_).empty());

    BOOST_CHECK(graph.getPredecessors(&endProcessor2->inport_).count(startProcessor)   &&
                !graph.getPredecessors(&endProcessor2->inport_).count(middleProcessor) &&
                !graph.getPredecessors(&endProcessor2->inport_).count(endProcessor2)   &&
                !graph.getPredecessors(&endProcessor2->inport_).count(endProcessor)
        );
    BOOST_CHECK(graph.getPredecessors(&endProcessor2->inport2_).empty());

    BOOST_CHECK(!graph.getSuccessors(&startProcessor->outport_).count(startProcessor)  &&
                graph.getSuccessors(&startProcessor->outport_).count(middleProcessor)  &&
                graph.getSuccessors(&startProcessor->outport_).count(endProcessor)     &&
                !graph.getSuccessors(&startProcessor->outport_).count(endProcessor2)
        );
    BOOST_CHECK(!graph.getSuccessors(&startProcessor->outport2_).count(startProcessor)  &&
                !graph.getSuccessors(&startProcessor->outport2_).count(middleProcessor) &&
                !graph.getSuccessors(&startProcessor->outport2_).count(endProcessor)    &&
                graph.getSuccessors(&startProcessor->outport2_).count(endProcessor2)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, &endProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, &endProcessor->inport2_));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, endProcessor));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport_, &endProcessor2->inport_));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport_, endProcessor2));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, endProcessor));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport2_, endProcessor2));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport2_, &endProcessor2->inport_));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, &endProcessor->inport_));

    BOOST_CHECK(graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport2_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport2_, &endProcessor2->inport_));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 4);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, startProcessor, endProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor));
}

// Network = (StartProcessor->MiddleProcessor->EndProcessor)
//            StartProcessor2-^
BOOST_AUTO_TEST_CASE(merge)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor);
    processors.push_back(startProcessor2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outport_, middleProcessor->inport2_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, endProcessor->inport_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(startProcessor2));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // check predecessors/successors
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)    &&
                graph.getPredecessors(endProcessor).count(startProcessor2)   &&
                graph.getPredecessors(endProcessor).count(middleProcessor)   &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)    &&
                graph.getSuccessors(startProcessor).count(middleProcessor)   &&
                graph.getSuccessors(startProcessor).count(endProcessor)      &&
                !graph.getSuccessors(startProcessor).count(startProcessor2)
        );
    BOOST_CHECK(graph.getSuccessors(startProcessor2).count(startProcessor2)  &&
                graph.getSuccessors(startProcessor2).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor2).count(endProcessor)     &&
                !graph.getSuccessors(startProcessor2).count(startProcessor)
        );
    std::set<Processor*> startProcessors;
    startProcessors.insert(startProcessor);
    startProcessors.insert(startProcessor2);
    BOOST_CHECK(graph.getSuccessors(startProcessors).count(startProcessor2)  &&
                graph.getSuccessors(startProcessors).count(middleProcessor)  &&
                graph.getSuccessors(startProcessors).count(endProcessor)     &&
                graph.getSuccessors(startProcessors).count(startProcessor)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor2, endProcessor));

    // port predecessors/successors
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport_).count(startProcessor)   &&
                graph.getPredecessors(&endProcessor->inport_).count(startProcessor2)  &&
                graph.getPredecessors(&endProcessor->inport_).count(middleProcessor)  &&
                !graph.getPredecessors(&endProcessor->inport_).count(endProcessor)
        );
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport2_).empty());
    BOOST_CHECK(graph.getPredecessors(&middleProcessor->inport_).count(startProcessor)    &&
                !graph.getPredecessors(&middleProcessor->inport_).count(startProcessor2)  &&
                !graph.getPredecessors(&middleProcessor->inport_).count(middleProcessor)  &&
                !graph.getPredecessors(&middleProcessor->inport_).count(endProcessor)
        );
    BOOST_CHECK(!graph.getPredecessors(&middleProcessor->inport2_).count(startProcessor)  &&
                graph.getPredecessors(&middleProcessor->inport2_).count(startProcessor2)  &&
                !graph.getPredecessors(&middleProcessor->inport2_).count(middleProcessor) &&
                !graph.getPredecessors(&middleProcessor->inport2_).count(endProcessor)
        );


    BOOST_CHECK(!graph.getSuccessors(&startProcessor->outport_).count(startProcessor)   &&
                graph.getSuccessors(&startProcessor->outport_).count(middleProcessor)   &&
                graph.getSuccessors(&startProcessor->outport_).count(endProcessor)      &&
                !graph.getSuccessors(&startProcessor->outport_).count(startProcessor2)
        );
    BOOST_CHECK(graph.getSuccessors(&startProcessor->outport2_).empty());
    BOOST_CHECK(!graph.getSuccessors(&startProcessor2->outport_).count(startProcessor2) &&
                graph.getSuccessors(&startProcessor2->outport_).count(middleProcessor)  &&
                graph.getSuccessors(&startProcessor2->outport_).count(endProcessor)     &&
                !graph.getSuccessors(&startProcessor2->outport_).count(startProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(&startProcessor2->outport2_).empty());

    BOOST_CHECK(graph.isSuccessor(startProcessor, &endProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, &endProcessor->inport2_));
    BOOST_CHECK(graph.isSuccessor(startProcessor, &endProcessor->inport_));
    BOOST_CHECK(graph.isSuccessor(startProcessor, &middleProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor2, &middleProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, &middleProcessor->inport2_));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, &middleProcessor->inport2_));

    BOOST_CHECK(graph.isSuccessor(&startProcessor2->outport_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor2->outport2_, endProcessor));

    BOOST_CHECK(graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport2_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport2_));
    BOOST_CHECK(graph.isPathElement(middleProcessor, &startProcessor2->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor2->outport2_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor2->outport_, &endProcessor->inport2_));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 4);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, startProcessor2, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor));
}

// Network = (StartProcessor->MiddleProcessor->EndProcessor )
//           (                              \->EndProcessor2)
BOOST_AUTO_TEST_CASE(split)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor2);
    processors.push_back(endProcessor);
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport2_, endProcessor2->inport_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor2));

    // check predecessors/successors
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(middleProcessor)  &&
                graph.getPredecessors(endProcessor).count(endProcessor)     &&
                !graph.getPredecessors(endProcessor).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(endProcessor2).count(startProcessor)  &&
                graph.getPredecessors(endProcessor2).count(middleProcessor) &&
                graph.getPredecessors(endProcessor2).count(endProcessor2)   &&
                !graph.getPredecessors(endProcessor2).count(endProcessor)
    );
    std::set<Processor*> endProcessors;
    endProcessors.insert(endProcessor);
    endProcessors.insert(endProcessor2);
    BOOST_CHECK(graph.getPredecessors(endProcessors).count(startProcessor)  &&
                graph.getPredecessors(endProcessors).count(middleProcessor) &&
                graph.getPredecessors(endProcessors).count(endProcessor)    &&
                graph.getPredecessors(endProcessors).count(endProcessor2)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor).count(endProcessor)     &&
                graph.getSuccessors(startProcessor).count(endProcessor2)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor2));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor2));

    // port predecessors/successors
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport_).count(startProcessor)   &&
                graph.getPredecessors(&endProcessor->inport_).count(middleProcessor)  &&
                !graph.getPredecessors(&endProcessor->inport_).count(endProcessor)    &&
                !graph.getPredecessors(&endProcessor->inport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport2_).empty());
    BOOST_CHECK(graph.getPredecessors(&endProcessor2->inport_).count(startProcessor)  &&
                graph.getPredecessors(&endProcessor2->inport_).count(middleProcessor) &&
                !graph.getPredecessors(&endProcessor2->inport_).count(endProcessor2)  &&
                !graph.getPredecessors(&endProcessor2->inport_).count(endProcessor)
        );
    BOOST_CHECK(graph.getPredecessors(&endProcessor2->inport2_).empty());

    BOOST_CHECK(!graph.getSuccessors(&startProcessor->outport_).count(startProcessor) &&
                graph.getSuccessors(&startProcessor->outport_).count(middleProcessor) &&
                graph.getSuccessors(&startProcessor->outport_).count(endProcessor)    &&
                graph.getSuccessors(&startProcessor->outport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getSuccessors(&startProcessor->outport2_).empty());
    BOOST_CHECK(!graph.getSuccessors(&middleProcessor->outport_).count(startProcessor)  &&
                !graph.getSuccessors(&middleProcessor->outport_).count(middleProcessor) &&
                graph.getSuccessors(&middleProcessor->outport_).count(endProcessor)     &&
                !graph.getSuccessors(&middleProcessor->outport_).count(endProcessor2)
        );
    BOOST_CHECK(!graph.getSuccessors(&middleProcessor->outport2_).count(startProcessor)  &&
                !graph.getSuccessors(&middleProcessor->outport2_).count(middleProcessor) &&
                !graph.getSuccessors(&middleProcessor->outport2_).count(endProcessor)    &&
                graph.getSuccessors(&middleProcessor->outport2_).count(endProcessor2)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, &endProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, &endProcessor->inport2_));
    BOOST_CHECK(graph.isSuccessor(startProcessor, &endProcessor2->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, &endProcessor2->inport2_));

    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, endProcessor));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, endProcessor2));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, endProcessor2));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, endProcessor2));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, endProcessor2));
    BOOST_CHECK(graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport2_));
    BOOST_CHECK(graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor2->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor2->inport2_));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 4);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor2));
}

// Network = (StartProcessor->MiddleProcessor->EndProcessor)
//           (             \->MiddleProcessor2-^           )
BOOST_AUTO_TEST_CASE(splitMerge)
{
    vector<Processor*> processors;
    processors.push_back(middleProcessor);
    processors.push_back(middleProcessor2);
    processors.push_back(endProcessor);
    processors.push_back(startProcessor);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(startProcessor->outport2_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, endProcessor->inport2_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor2));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // check predecessors/successors
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(middleProcessor)  &&
                graph.getPredecessors(endProcessor).count(middleProcessor2) &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor).count(middleProcessor2) &&
                graph.getSuccessors(startProcessor).count(endProcessor)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor, endProcessor));

    // port predecessors/successors
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport_).count(startProcessor)    &&
                graph.getPredecessors(&endProcessor->inport_).count(middleProcessor)   &&
                !graph.getPredecessors(&endProcessor->inport_).count(middleProcessor2) &&
                !graph.getPredecessors(&endProcessor->inport_).count(endProcessor)
        );
    BOOST_CHECK(graph.getPredecessors(&endProcessor->inport2_).count(startProcessor)   &&
                !graph.getPredecessors(&endProcessor->inport2_).count(middleProcessor) &&
                graph.getPredecessors(&endProcessor->inport2_).count(middleProcessor2) &&
                !graph.getPredecessors(&endProcessor->inport2_).count(endProcessor)
        );

    BOOST_CHECK(!graph.getSuccessors(&startProcessor->outport_).count(startProcessor)   &&
                graph.getSuccessors(&startProcessor->outport_).count(middleProcessor)   &&
                !graph.getSuccessors(&startProcessor->outport_).count(middleProcessor2) &&
                graph.getSuccessors(&startProcessor->outport_).count(endProcessor)
        );
    BOOST_CHECK(!graph.getSuccessors(&startProcessor->outport2_).count(startProcessor)  &&
                !graph.getSuccessors(&startProcessor->outport2_).count(middleProcessor) &&
                graph.getSuccessors(&startProcessor->outport2_).count(middleProcessor2) &&
                graph.getSuccessors(&startProcessor->outport2_).count(endProcessor)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, &endProcessor->inport_));
    BOOST_CHECK(graph.isSuccessor(startProcessor, &endProcessor->inport2_));
    BOOST_CHECK(graph.isSuccessor(middleProcessor, &endProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(middleProcessor, &endProcessor->inport2_));
    BOOST_CHECK(!graph.isSuccessor(middleProcessor2, &endProcessor->inport_));
    BOOST_CHECK(graph.isSuccessor(middleProcessor2, &endProcessor->inport2_));

    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, endProcessor));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport2_, endProcessor));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport_, middleProcessor));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport2_, middleProcessor));
    BOOST_CHECK(!graph.isSuccessor(&startProcessor->outport_, middleProcessor2));
    BOOST_CHECK(graph.isSuccessor(&startProcessor->outport2_, middleProcessor2));

    BOOST_CHECK(graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport2_, &endProcessor->inport2_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, &startProcessor->outport2_, &endProcessor->inport_));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, &startProcessor->outport2_, &endProcessor->inport2_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, &startProcessor->outport_, &endProcessor->inport2_));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 4);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor2, endProcessor));
}

// Network = (StartProcessor ->MiddleProcessor ->EndProcessor )
//            StartProcessor2->MiddleProcessor2->EndProcessor2)
BOOST_AUTO_TEST_CASE(separatePipelines)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor2);
    processors.push_back(middleProcessor);
    processors.push_back(endProcessor2);
    processors.push_back(startProcessor2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, endProcessor2->inport_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 6);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(startProcessor2));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor2));
    BOOST_CHECK(graph.containsProcessor(endProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor2));

    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)     &&
                graph.getPredecessors(endProcessor).count(middleProcessor)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)       &&
                !graph.getPredecessors(endProcessor).count(startProcessor2)   &&
                !graph.getPredecessors(endProcessor).count(middleProcessor2)  &&
                !graph.getPredecessors(endProcessor).count(endProcessor2)
        );
    BOOST_CHECK(!graph.getPredecessors(endProcessor2).count(startProcessor)   &&
                !graph.getPredecessors(endProcessor2).count(middleProcessor)  &&
                !graph.getPredecessors(endProcessor2).count(endProcessor)     &&
                graph.getPredecessors(endProcessor2).count(startProcessor2)   &&
                graph.getPredecessors(endProcessor2).count(middleProcessor2)  &&
                graph.getPredecessors(endProcessor2).count(endProcessor2)
        );
    std::set<Processor*> endProcessors;
    endProcessors.insert(endProcessor);
    endProcessors.insert(endProcessor2);
    BOOST_CHECK(graph.getPredecessors(endProcessors).count(startProcessor)    &&
                graph.getPredecessors(endProcessors).count(middleProcessor)   &&
                graph.getPredecessors(endProcessors).count(endProcessor)      &&
                graph.getPredecessors(endProcessors).count(startProcessor2)   &&
                graph.getPredecessors(endProcessors).count(middleProcessor2)  &&
                graph.getPredecessors(endProcessors).count(endProcessor2)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)     &&
                graph.getSuccessors(startProcessor).count(middleProcessor)    &&
                graph.getSuccessors(startProcessor).count(endProcessor)       &&
                !graph.getSuccessors(startProcessor).count(startProcessor2)   &&
                !graph.getSuccessors(startProcessor).count(middleProcessor2)  &&
                !graph.getSuccessors(startProcessor).count(endProcessor2)
        );
    BOOST_CHECK(!graph.getSuccessors(startProcessor2).count(startProcessor)   &&
                !graph.getSuccessors(startProcessor2).count(middleProcessor)  &&
                !graph.getSuccessors(startProcessor2).count(endProcessor)     &&
                graph.getSuccessors(startProcessor2).count(startProcessor2)   &&
                graph.getSuccessors(startProcessor2).count(middleProcessor2)  &&
                graph.getSuccessors(startProcessor2).count(endProcessor2)
        );
    std::set<Processor*> startProcessors;
    startProcessors.insert(startProcessor);
    startProcessors.insert(startProcessor2);
    BOOST_CHECK(graph.getSuccessors(startProcessors).count(startProcessor)    &&
                graph.getSuccessors(startProcessors).count(middleProcessor)   &&
                graph.getSuccessors(startProcessors).count(endProcessor)      &&
                graph.getSuccessors(startProcessors).count(startProcessor2)   &&
                graph.getSuccessors(startProcessors).count(middleProcessor2)  &&
                graph.getSuccessors(startProcessors).count(endProcessor2)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, endProcessor2));
    BOOST_CHECK(!graph.isSuccessor(startProcessor2, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, endProcessor2));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor2, endProcessor2));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, startProcessor, endProcessor));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor2, endProcessor2));
    BOOST_CHECK(graph.isPathElement(middleProcessor, &startProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, &startProcessor2->outport_, &endProcessor2->inport_));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 6);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor));
    BOOST_CHECK(beforeInVector(sorted, startProcessor2, middleProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor2, endProcessor2));
}

// Network = (StartProcessor ->MiddleProcessor ->EndProcessor )
//                                   v
//            StartProcessor2->MiddleProcessor2->EndProcessor2)
BOOST_AUTO_TEST_CASE(connectedPipelines)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor2);
    processors.push_back(middleProcessor);
    processors.push_back(endProcessor2);
    processors.push_back(startProcessor2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport2_, middleProcessor2->inport2_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, endProcessor2->inport_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 6);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(startProcessor2));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor2));
    BOOST_CHECK(graph.containsProcessor(endProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor2));

    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)     &&
                graph.getPredecessors(endProcessor).count(middleProcessor)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)       &&
                !graph.getPredecessors(endProcessor).count(startProcessor2)   &&
                !graph.getPredecessors(endProcessor).count(middleProcessor2)  &&
                !graph.getPredecessors(endProcessor).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(endProcessor2).count(startProcessor)    &&
                graph.getPredecessors(endProcessor2).count(middleProcessor)   &&
                !graph.getPredecessors(endProcessor2).count(endProcessor)     &&
                graph.getPredecessors(endProcessor2).count(startProcessor2)   &&
                graph.getPredecessors(endProcessor2).count(middleProcessor2)  &&
                graph.getPredecessors(endProcessor2).count(endProcessor2)
        );
    std::set<Processor*> endProcessors;
    endProcessors.insert(endProcessor);
    endProcessors.insert(endProcessor2);
    BOOST_CHECK(graph.getPredecessors(endProcessors).count(startProcessor)    &&
                graph.getPredecessors(endProcessors).count(middleProcessor)   &&
                graph.getPredecessors(endProcessors).count(endProcessor)      &&
                graph.getPredecessors(endProcessors).count(startProcessor2)   &&
                graph.getPredecessors(endProcessors).count(middleProcessor2)  &&
                graph.getPredecessors(endProcessors).count(endProcessor2)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)     &&
                graph.getSuccessors(startProcessor).count(middleProcessor)    &&
                graph.getSuccessors(startProcessor).count(endProcessor)       &&
                !graph.getSuccessors(startProcessor).count(startProcessor2)   &&
                graph.getSuccessors(startProcessor).count(middleProcessor2)   &&
                graph.getSuccessors(startProcessor).count(endProcessor2)
        );
    BOOST_CHECK(!graph.getSuccessors(startProcessor2).count(startProcessor)   &&
                !graph.getSuccessors(startProcessor2).count(middleProcessor)  &&
                !graph.getSuccessors(startProcessor2).count(endProcessor)     &&
                graph.getSuccessors(startProcessor2).count(startProcessor2)   &&
                graph.getSuccessors(startProcessor2).count(middleProcessor2)  &&
                graph.getSuccessors(startProcessor2).count(endProcessor2)
        );
    std::set<Processor*> startProcessors;
    startProcessors.insert(startProcessor);
    startProcessors.insert(startProcessor2);
    BOOST_CHECK(graph.getSuccessors(startProcessors).count(startProcessor)    &&
                graph.getSuccessors(startProcessors).count(middleProcessor)   &&
                graph.getSuccessors(startProcessors).count(endProcessor)      &&
                graph.getSuccessors(startProcessors).count(startProcessor2)   &&
                graph.getSuccessors(startProcessors).count(middleProcessor2)  &&
                graph.getSuccessors(startProcessors).count(endProcessor2)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor2));
    BOOST_CHECK(!graph.isSuccessor(startProcessor2, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, endProcessor2));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor2, endProcessor2));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor, endProcessor2));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor2, endProcessor2));

    // port predecessors/successors
    BOOST_CHECK(graph.getPredecessors(&middleProcessor->inport_).count(startProcessor)     &&
                !graph.getPredecessors(&middleProcessor->inport_).count(middleProcessor)   &&
                !graph.getPredecessors(&middleProcessor->inport_).count(endProcessor)      &&
                !graph.getPredecessors(&middleProcessor->inport_).count(startProcessor2)   &&
                !graph.getPredecessors(&middleProcessor->inport_).count(middleProcessor2)  &&
                !graph.getPredecessors(&middleProcessor->inport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(&middleProcessor->inport2_).empty());

    BOOST_CHECK(!graph.getPredecessors(&middleProcessor2->inport_).count(startProcessor)   &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(middleProcessor)  &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(endProcessor)     &&
                graph.getPredecessors(&middleProcessor2->inport_).count(startProcessor2)   &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(middleProcessor2) &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(&middleProcessor2->inport2_).count(startProcessor)    &&
                graph.getPredecessors(&middleProcessor2->inport2_).count(middleProcessor)   &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(endProcessor)     &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(startProcessor2)  &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(middleProcessor2) &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(endProcessor2)
        );

    BOOST_CHECK(!graph.getSuccessors(&middleProcessor->outport_).count(startProcessor)    &&
                !graph.getSuccessors(&middleProcessor->outport_).count(middleProcessor)   &&
                graph.getSuccessors(&middleProcessor->outport_).count(endProcessor)       &&
                !graph.getSuccessors(&middleProcessor->outport_).count(startProcessor2)   &&
                !graph.getSuccessors(&middleProcessor->outport_).count(middleProcessor2)  &&
                !graph.getSuccessors(&middleProcessor->outport_).count(endProcessor2)
        );
    BOOST_CHECK(!graph.getSuccessors(&middleProcessor->outport2_).count(startProcessor)   &&
                !graph.getSuccessors(&middleProcessor->outport2_).count(middleProcessor)  &&
                !graph.getSuccessors(&middleProcessor->outport2_).count(endProcessor)     &&
                !graph.getSuccessors(&middleProcessor->outport2_).count(startProcessor2)  &&
                graph.getSuccessors(&middleProcessor->outport2_).count(middleProcessor2)  &&
                graph.getSuccessors(&middleProcessor->outport2_).count(endProcessor2)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, &middleProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor2, &middleProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, &middleProcessor->inport2_));

    BOOST_CHECK(!graph.isSuccessor(startProcessor, &middleProcessor2->inport_));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, &middleProcessor2->inport_));
    BOOST_CHECK(graph.isSuccessor(startProcessor, &middleProcessor2->inport2_));
    BOOST_CHECK(graph.isSuccessor(middleProcessor, &middleProcessor2->inport2_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor2, &middleProcessor2->inport2_));

    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, endProcessor2));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport2_, endProcessor));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport2_, endProcessor2));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, middleProcessor2));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport2_, middleProcessor2));

    //BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, &middleProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, &middleProcessor->outport_, &endProcessor2->inport_));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, &middleProcessor->outport2_, &endProcessor2->inport_));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 6);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor));
    BOOST_CHECK(beforeInVector(sorted, startProcessor2, middleProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor2, endProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, middleProcessor2));
}

// Network = (StartProcessor ->MiddleProcessor ->EndProcessor )
//                                   v    /-------^
//            StartProcessor2->MiddleProcessor2->EndProcessor2)
BOOST_AUTO_TEST_CASE(crossconnectedPipelines)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor2);
    processors.push_back(middleProcessor);
    processors.push_back(endProcessor2);
    processors.push_back(startProcessor2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport2_, middleProcessor2->inport2_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, endProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport2_, endProcessor->inport2_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 6);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(startProcessor2));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor2));
    BOOST_CHECK(graph.containsProcessor(endProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor2));

    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)     &&
                graph.getPredecessors(endProcessor).count(middleProcessor)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)       &&
                graph.getPredecessors(endProcessor).count(startProcessor2)    &&
                graph.getPredecessors(endProcessor).count(middleProcessor2)   &&
                !graph.getPredecessors(endProcessor).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(endProcessor2).count(startProcessor)    &&
                graph.getPredecessors(endProcessor2).count(middleProcessor)   &&
                !graph.getPredecessors(endProcessor2).count(endProcessor)     &&
                graph.getPredecessors(endProcessor2).count(startProcessor2)   &&
                graph.getPredecessors(endProcessor2).count(middleProcessor2)  &&
                graph.getPredecessors(endProcessor2).count(endProcessor2)
        );
    std::set<Processor*> endProcessors;
    endProcessors.insert(endProcessor);
    endProcessors.insert(endProcessor2);
    BOOST_CHECK(graph.getPredecessors(endProcessors).count(startProcessor)    &&
                graph.getPredecessors(endProcessors).count(middleProcessor)   &&
                graph.getPredecessors(endProcessors).count(endProcessor)      &&
                graph.getPredecessors(endProcessors).count(startProcessor2)   &&
                graph.getPredecessors(endProcessors).count(middleProcessor2)  &&
                graph.getPredecessors(endProcessors).count(endProcessor2)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)     &&
                graph.getSuccessors(startProcessor).count(middleProcessor)    &&
                graph.getSuccessors(startProcessor).count(endProcessor)       &&
                !graph.getSuccessors(startProcessor).count(startProcessor2)   &&
                graph.getSuccessors(startProcessor).count(middleProcessor2)   &&
                graph.getSuccessors(startProcessor).count(endProcessor2)
        );
    BOOST_CHECK(!graph.getSuccessors(startProcessor2).count(startProcessor)   &&
                !graph.getSuccessors(startProcessor2).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor2).count(endProcessor)      &&
                graph.getSuccessors(startProcessor2).count(startProcessor2)   &&
                graph.getSuccessors(startProcessor2).count(middleProcessor2)  &&
                graph.getSuccessors(startProcessor2).count(endProcessor2)
        );
    std::set<Processor*> startProcessors;
    startProcessors.insert(startProcessor);
    startProcessors.insert(startProcessor2);
    BOOST_CHECK(graph.getSuccessors(startProcessors).count(startProcessor)    &&
                graph.getSuccessors(startProcessors).count(middleProcessor)   &&
                graph.getSuccessors(startProcessors).count(endProcessor)      &&
                graph.getSuccessors(startProcessors).count(startProcessor2)   &&
                graph.getSuccessors(startProcessors).count(middleProcessor2)  &&
                graph.getSuccessors(startProcessors).count(endProcessor2)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor2));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, endProcessor2));

    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor2));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor, endProcessor2));

    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor2, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor2, endProcessor));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor2, endProcessor2));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor2, endProcessor2));

    // port predecessors/successors
    BOOST_CHECK(graph.getPredecessors(&middleProcessor->inport_).count(startProcessor)     &&
                !graph.getPredecessors(&middleProcessor->inport_).count(middleProcessor)   &&
                !graph.getPredecessors(&middleProcessor->inport_).count(endProcessor)      &&
                !graph.getPredecessors(&middleProcessor->inport_).count(startProcessor2)   &&
                !graph.getPredecessors(&middleProcessor->inport_).count(middleProcessor2)  &&
                !graph.getPredecessors(&middleProcessor->inport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(&middleProcessor->inport2_).empty());

    BOOST_CHECK(!graph.getPredecessors(&middleProcessor2->inport_).count(startProcessor)   &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(middleProcessor)  &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(endProcessor)     &&
                graph.getPredecessors(&middleProcessor2->inport_).count(startProcessor2)   &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(middleProcessor2) &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(&middleProcessor2->inport2_).count(startProcessor)    &&
                graph.getPredecessors(&middleProcessor2->inport2_).count(middleProcessor)   &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(endProcessor)     &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(startProcessor2)  &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(middleProcessor2) &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(endProcessor2)
        );

    BOOST_CHECK(!graph.getSuccessors(&middleProcessor->outport_).count(startProcessor)    &&
                !graph.getSuccessors(&middleProcessor->outport_).count(middleProcessor)   &&
                graph.getSuccessors(&middleProcessor->outport_).count(endProcessor)       &&
                !graph.getSuccessors(&middleProcessor->outport_).count(startProcessor2)   &&
                !graph.getSuccessors(&middleProcessor->outport_).count(middleProcessor2)  &&
                !graph.getSuccessors(&middleProcessor->outport_).count(endProcessor2)
        );
    BOOST_CHECK(!graph.getSuccessors(&middleProcessor->outport2_).count(startProcessor)   &&
                !graph.getSuccessors(&middleProcessor->outport2_).count(middleProcessor)  &&
                graph.getSuccessors(&middleProcessor->outport2_).count(endProcessor)     &&
                !graph.getSuccessors(&middleProcessor->outport2_).count(startProcessor2)  &&
                graph.getSuccessors(&middleProcessor->outport2_).count(middleProcessor2)  &&
                graph.getSuccessors(&middleProcessor->outport2_).count(endProcessor2)
        );


    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, endProcessor2));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport2_, endProcessor));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport2_, endProcessor2));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, middleProcessor2));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport2_, middleProcessor2));

    //BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    //BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor2, endProcessor2));
    //BOOST_CHECK(!graph.isPathElement(middleProcessor2, startProcessor, endProcessor));
    //BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor, endProcessor2));
    //BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor2, endProcessor2));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 6);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor));
    BOOST_CHECK(beforeInVector(sorted, startProcessor2, middleProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor2, endProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, middleProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor2, endProcessor2));
}

// Network = (StartProcessor ->MiddleProcessor ->EndProcessor )
//                                   v    /-------^
//            StartProcessor2->MiddleProcessor2->EndProcessor2)
BOOST_AUTO_TEST_CASE(crossconnectedPipelinesMultiConnections)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(startProcessor);
    processors.push_back(middleProcessor2);
    processors.push_back(middleProcessor);
    processors.push_back(endProcessor2);
    processors.push_back(startProcessor2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, middleProcessor2->inport2_));  //<<<!!!
    BOOST_CHECK(network->connectPorts(startProcessor2->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, endProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, endProcessor->inport2_));     //<<<!!!

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 6);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(startProcessor2));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor2));
    BOOST_CHECK(graph.containsProcessor(endProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor2));

    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)     &&
                graph.getPredecessors(endProcessor).count(middleProcessor)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)       &&
                graph.getPredecessors(endProcessor).count(startProcessor2)    &&
                graph.getPredecessors(endProcessor).count(middleProcessor2)   &&
                !graph.getPredecessors(endProcessor).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(endProcessor2).count(startProcessor)    &&
                graph.getPredecessors(endProcessor2).count(middleProcessor)   &&
                !graph.getPredecessors(endProcessor2).count(endProcessor)     &&
                graph.getPredecessors(endProcessor2).count(startProcessor2)   &&
                graph.getPredecessors(endProcessor2).count(middleProcessor2)  &&
                graph.getPredecessors(endProcessor2).count(endProcessor2)
        );
    std::set<Processor*> endProcessors;
    endProcessors.insert(endProcessor);
    endProcessors.insert(endProcessor2);
    BOOST_CHECK(graph.getPredecessors(endProcessors).count(startProcessor)    &&
                graph.getPredecessors(endProcessors).count(middleProcessor)   &&
                graph.getPredecessors(endProcessors).count(endProcessor)      &&
                graph.getPredecessors(endProcessors).count(startProcessor2)   &&
                graph.getPredecessors(endProcessors).count(middleProcessor2)  &&
                graph.getPredecessors(endProcessors).count(endProcessor2)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)     &&
                graph.getSuccessors(startProcessor).count(middleProcessor)    &&
                graph.getSuccessors(startProcessor).count(endProcessor)       &&
                !graph.getSuccessors(startProcessor).count(startProcessor2)   &&
                graph.getSuccessors(startProcessor).count(middleProcessor2)   &&
                graph.getSuccessors(startProcessor).count(endProcessor2)
        );
    BOOST_CHECK(!graph.getSuccessors(startProcessor2).count(startProcessor)   &&
                !graph.getSuccessors(startProcessor2).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor2).count(endProcessor)      &&
                graph.getSuccessors(startProcessor2).count(startProcessor2)   &&
                graph.getSuccessors(startProcessor2).count(middleProcessor2)  &&
                graph.getSuccessors(startProcessor2).count(endProcessor2)
        );
    std::set<Processor*> startProcessors;
    startProcessors.insert(startProcessor);
    startProcessors.insert(startProcessor2);
    BOOST_CHECK(graph.getSuccessors(startProcessors).count(startProcessor)    &&
                graph.getSuccessors(startProcessors).count(middleProcessor)   &&
                graph.getSuccessors(startProcessors).count(endProcessor)      &&
                graph.getSuccessors(startProcessors).count(startProcessor2)   &&
                graph.getSuccessors(startProcessors).count(middleProcessor2)  &&
                graph.getSuccessors(startProcessors).count(endProcessor2)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor2));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, endProcessor2));

    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor2));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor, endProcessor2));

    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor2, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor2, endProcessor));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, startProcessor2, endProcessor2));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, startProcessor2, endProcessor2));

    // port predecessors/successors
    BOOST_CHECK(graph.getPredecessors(&middleProcessor->inport_).count(startProcessor)     &&
                !graph.getPredecessors(&middleProcessor->inport_).count(middleProcessor)   &&
                !graph.getPredecessors(&middleProcessor->inport_).count(endProcessor)      &&
                !graph.getPredecessors(&middleProcessor->inport_).count(startProcessor2)   &&
                !graph.getPredecessors(&middleProcessor->inport_).count(middleProcessor2)  &&
                !graph.getPredecessors(&middleProcessor->inport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(&middleProcessor->inport2_).empty());

    BOOST_CHECK(!graph.getPredecessors(&middleProcessor2->inport_).count(startProcessor)   &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(middleProcessor)  &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(endProcessor)     &&
                graph.getPredecessors(&middleProcessor2->inport_).count(startProcessor2)   &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(middleProcessor2) &&
                !graph.getPredecessors(&middleProcessor2->inport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(&middleProcessor2->inport2_).count(startProcessor)    &&
                graph.getPredecessors(&middleProcessor2->inport2_).count(middleProcessor)   &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(endProcessor)     &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(startProcessor2)  &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(middleProcessor2) &&
                !graph.getPredecessors(&middleProcessor2->inport2_).count(endProcessor2)
        );

    BOOST_CHECK(!graph.getSuccessors(&middleProcessor->outport_).count(startProcessor)    &&
                !graph.getSuccessors(&middleProcessor->outport_).count(middleProcessor)   &&
                graph.getSuccessors(&middleProcessor->outport_).count(endProcessor)       &&
                !graph.getSuccessors(&middleProcessor->outport_).count(startProcessor2)   &&
                graph.getSuccessors(&middleProcessor->outport_).count(middleProcessor2)   &&
                graph.getSuccessors(&middleProcessor->outport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getSuccessors(&middleProcessor->outport2_).empty());

    BOOST_CHECK(!graph.getSuccessors(&middleProcessor2->outport_).count(startProcessor)    &&
                !graph.getSuccessors(&middleProcessor2->outport_).count(middleProcessor)   &&
                graph.getSuccessors(&middleProcessor2->outport_).count(endProcessor)       &&
                !graph.getSuccessors(&middleProcessor2->outport_).count(startProcessor2)   &&
                !graph.getSuccessors(&middleProcessor2->outport_).count(middleProcessor2)  &&
                graph.getSuccessors(&middleProcessor2->outport_).count(endProcessor2)
        );
    BOOST_CHECK(graph.getSuccessors(&middleProcessor2->outport2_).empty());

    BOOST_CHECK(graph.isSuccessor(startProcessor, &middleProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor2, &middleProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor, &middleProcessor->inport2_));

    BOOST_CHECK(!graph.isSuccessor(startProcessor, &middleProcessor2->inport_));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, &middleProcessor2->inport_));
    BOOST_CHECK(graph.isSuccessor(startProcessor, &middleProcessor2->inport2_));
    BOOST_CHECK(graph.isSuccessor(middleProcessor, &middleProcessor2->inport2_));
    BOOST_CHECK(!graph.isSuccessor(startProcessor2, &middleProcessor2->inport2_));

    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, endProcessor));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, endProcessor2));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport2_, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport2_, endProcessor2));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, middleProcessor2));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport2_, middleProcessor2));

    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, &endProcessor->inport2_));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, &endProcessor2->inport_));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, &endProcessor2->inport2_));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, &middleProcessor->inport_));
    BOOST_CHECK(!graph.isSuccessor(&middleProcessor->outport_, &middleProcessor2->inport_));
    BOOST_CHECK(graph.isSuccessor(&middleProcessor->outport_, &middleProcessor2->inport2_));

    BOOST_CHECK(graph.isPathElement(middleProcessor2, &middleProcessor->outport_, &endProcessor->inport2_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, &middleProcessor->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, &middleProcessor->outport2_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, &middleProcessor2->outport_, &endProcessor->inport_));
    BOOST_CHECK(!graph.isPathElement(endProcessor, &middleProcessor->outport_, &endProcessor->inport_));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 6);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor));
    BOOST_CHECK(beforeInVector(sorted, startProcessor2, middleProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor2, endProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, middleProcessor2));
    BOOST_CHECK(beforeInVector(sorted, middleProcessor2, endProcessor2));
}

// Network = (StartProcessor->(LoopInitiator->MiddleProcessor->LoopFinalizer)x3->EndProcessor)
BOOST_AUTO_TEST_CASE(simpleLoop)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(loopInitiator);
    processors.push_back(middleProcessor);
    processors.push_back(startProcessor);
    processors.push_back(loopFinalizer);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutport_, loopInitiator->loopInport_));
    loopInitiator->loopInport_.setNumLoopIterations(3);

    NetworkGraph graph(processors);

    //BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(loopInitiator));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // processors are connected => processors are predecessors/successors of each other
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(loopInitiator)    &&
                graph.getPredecessors(endProcessor).count(middleProcessor)  &&
                graph.getPredecessors(endProcessor).count(loopFinalizer)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(loopInitiator)    &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor).count(loopFinalizer)    &&
                graph.getSuccessors(startProcessor).count(endProcessor)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(loopInitiator, loopFinalizer));
    BOOST_CHECK(graph.isSuccessor(loopInitiator, middleProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, loopInitiator, loopFinalizer));
    BOOST_CHECK(graph.isSuccessor(middleProcessor, loopInitiator));                  //< due to loop unrolling
    BOOST_CHECK(graph.isPathElement(middleProcessor, loopFinalizer, loopInitiator)); //< due to loop unrolling

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 11);
    BOOST_CHECK(sorted.at(0) == startProcessor);
    BOOST_CHECK(sorted.at(1) == loopInitiator);
    BOOST_CHECK(sorted.at(2) == middleProcessor);
    BOOST_CHECK(sorted.at(3) == loopFinalizer);
    BOOST_CHECK(sorted.at(4) == loopInitiator);
    BOOST_CHECK(sorted.at(5) == middleProcessor);
    BOOST_CHECK(sorted.at(6) == loopFinalizer);
    BOOST_CHECK(sorted.at(7) == loopInitiator);
    BOOST_CHECK(sorted.at(8) == middleProcessor);
    BOOST_CHECK(sorted.at(9) == loopFinalizer);
    BOOST_CHECK(sorted.at(10) == endProcessor);
}

// Network = (StartProcessor->(LoopInitiator->MiddleProcessor LoopFinalizer)x3->EndProcessor)
BOOST_AUTO_TEST_CASE(unconnectedLoop)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(loopInitiator);
    processors.push_back(middleProcessor);
    processors.push_back(startProcessor);
    processors.push_back(loopFinalizer);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, middleProcessor->inport_));
    //BOOST_CHECK(network->connectPorts(middleProcessor->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutport_, loopInitiator->loopInport_));
    loopInitiator->loopInport_.setNumLoopIterations(3);

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 5);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(loopInitiator));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // check predecessors/successors
    BOOST_CHECK(!graph.getPredecessors(endProcessor).count(startProcessor)  &&
                !graph.getPredecessors(endProcessor).count(loopInitiator)   &&
                !graph.getPredecessors(endProcessor).count(middleProcessor) &&
                graph.getPredecessors(endProcessor).count(loopFinalizer)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(loopInitiator)    &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                !graph.getSuccessors(startProcessor).count(loopFinalizer)   &&
                !graph.getSuccessors(startProcessor).count(endProcessor)
        );

    BOOST_CHECK(!graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(!graph.isSuccessor(loopInitiator, loopFinalizer));
    BOOST_CHECK(graph.isSuccessor(loopInitiator, middleProcessor));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, loopInitiator, loopFinalizer));
    BOOST_CHECK(!graph.isSuccessor(middleProcessor, loopInitiator));                  //< loop not unrolled
    BOOST_CHECK(!graph.isPathElement(middleProcessor, loopFinalizer, loopInitiator)); //< loop not unrolled

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 5);
    BOOST_CHECK(beforeInVector(sorted, startProcessor, loopInitiator));
    BOOST_CHECK(beforeInVector(sorted, loopInitiator, middleProcessor));
    BOOST_CHECK(beforeInVector(sorted, loopFinalizer, endProcessor));
}

// Network = (StartProcessor->LoopInitiator->MiddleProcessor->LoopFinalizer->EndProcessor)
BOOST_AUTO_TEST_CASE(inactiveLoop)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(loopInitiator);
    processors.push_back(middleProcessor);
    processors.push_back(startProcessor);
    processors.push_back(loopFinalizer);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));

    NetworkGraph graph(processors);

    BOOST_CHECK(graph.getProcessors().size() == 5);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(loopInitiator));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // processors are connected => processors are predecessors/successors of each other
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(loopInitiator)    &&
                graph.getPredecessors(endProcessor).count(middleProcessor)  &&
                graph.getPredecessors(endProcessor).count(loopFinalizer)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(loopInitiator)    &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor).count(loopFinalizer)    &&
                graph.getSuccessors(startProcessor).count(endProcessor)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(loopInitiator, loopFinalizer));
    BOOST_CHECK(graph.isSuccessor(loopInitiator, middleProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, loopInitiator, loopFinalizer));
    BOOST_CHECK(!graph.isSuccessor(middleProcessor, loopInitiator));                    //< loop not active
    BOOST_CHECK(!graph.isPathElement(middleProcessor, loopFinalizer, loopInitiator));   //< loop not active

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 5);
    BOOST_CHECK(sorted.at(0) == startProcessor);
    BOOST_CHECK(sorted.at(1) == loopInitiator);
    BOOST_CHECK(sorted.at(2) == middleProcessor);
    BOOST_CHECK(sorted.at(3) == loopFinalizer);
    BOOST_CHECK(sorted.at(4) == endProcessor);
}

// Network = (StartProcessor->(LoopInitiator->LoopFinalizer)x3->EndProcessor)
BOOST_AUTO_TEST_CASE(emptyLoop)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(loopInitiator);
    processors.push_back(startProcessor);
    processors.push_back(loopFinalizer);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutport_, loopInitiator->loopInport_));
    loopInitiator->loopInport_.setNumLoopIterations(3);

    NetworkGraph graph(processors);

    //BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(loopInitiator));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // processors are connected => processors are predecessors/successors of each other
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(loopInitiator)    &&
                graph.getPredecessors(endProcessor).count(loopFinalizer)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(loopInitiator)    &&
                graph.getSuccessors(startProcessor).count(loopFinalizer)    &&
                graph.getSuccessors(startProcessor).count(endProcessor)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(loopInitiator, loopFinalizer));
    BOOST_CHECK(graph.isSuccessor(loopFinalizer, loopInitiator));   //< loop connection

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 8);
    BOOST_CHECK(sorted.at(0) == startProcessor);
    BOOST_CHECK(sorted.at(1) == loopInitiator);
    BOOST_CHECK(sorted.at(2) == loopFinalizer);
    BOOST_CHECK(sorted.at(3) == loopInitiator);
    BOOST_CHECK(sorted.at(4) == loopFinalizer);
    BOOST_CHECK(sorted.at(5) == loopInitiator);
    BOOST_CHECK(sorted.at(6) == loopFinalizer);
    BOOST_CHECK(sorted.at(7) == endProcessor);
}

// Network = (StartProcessor->(LoopInitiator->(LoopInitiator2->MiddleProcessor2->LoopFinalizer2)x3->MiddleProcessor->LoopFinalizer)x2->EndProcessor)
BOOST_AUTO_TEST_CASE(nestedLoops)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(loopInitiator);
    processors.push_back(loopInitiator2);
    processors.push_back(middleProcessor);
    processors.push_back(middleProcessor2);
    processors.push_back(startProcessor);
    processors.push_back(loopFinalizer);
    processors.push_back(loopFinalizer2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, loopInitiator2->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator2->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, loopFinalizer2->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer2->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));

    BOOST_CHECK(network->connectPorts(loopFinalizer2->loopOutport_, loopInitiator2->loopInport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutport_, loopInitiator->loopInport_));

    loopInitiator2->loopInport_.setNumLoopIterations(3);
    loopInitiator->loopInport_.setNumLoopIterations(2);

    NetworkGraph graph(processors);

    //BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(loopInitiator));
    BOOST_CHECK(graph.containsProcessor(loopInitiator2));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor2));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer2));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // check predecessors/successors
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(loopInitiator)    &&
                graph.getPredecessors(endProcessor).count(loopInitiator2)   &&
                graph.getPredecessors(endProcessor).count(middleProcessor)  &&
                graph.getPredecessors(endProcessor).count(middleProcessor2) &&
                graph.getPredecessors(endProcessor).count(loopFinalizer)    &&
                graph.getPredecessors(endProcessor).count(loopFinalizer2)   &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(loopInitiator)    &&
                graph.getSuccessors(startProcessor).count(loopInitiator2)   &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor).count(middleProcessor2) &&
                graph.getSuccessors(startProcessor).count(loopFinalizer)    &&
                graph.getSuccessors(startProcessor).count(loopFinalizer2)   &&
                graph.getSuccessors(startProcessor).count(endProcessor)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    // loop features: due to the unrooling every processor in the loops is predecessor/successor of each other
    BOOST_CHECK(graph.isSuccessor(loopInitiator, loopFinalizer));
    BOOST_CHECK(graph.isSuccessor(loopInitiator, loopInitiator2));
    BOOST_CHECK(graph.isSuccessor(loopInitiator2, loopInitiator));
    BOOST_CHECK(graph.isSuccessor(loopFinalizer, loopInitiator));
    BOOST_CHECK(graph.isSuccessor(loopFinalizer2, loopInitiator2));
    BOOST_CHECK(graph.isSuccessor(loopFinalizer2, loopInitiator));
    BOOST_CHECK(graph.isPathElement(middleProcessor, loopInitiator, loopFinalizer));
    BOOST_CHECK(graph.isPathElement(middleProcessor, loopInitiator, loopInitiator2));
    BOOST_CHECK(graph.isPathElement(middleProcessor, loopInitiator2, loopInitiator));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, loopFinalizer, loopInitiator));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, loopFinalizer2, loopInitiator2));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, loopFinalizer2, loopInitiator));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 26);
    BOOST_CHECK(sorted.at(0) == startProcessor);
        // outer loop
        BOOST_CHECK(sorted.at(1) == loopInitiator);
            // inner loop
            BOOST_CHECK(sorted.at(2) == loopInitiator2);
            BOOST_CHECK(sorted.at(3) == middleProcessor2);
            BOOST_CHECK(sorted.at(4) == loopFinalizer2);
            BOOST_CHECK(sorted.at(5) == loopInitiator2);
            BOOST_CHECK(sorted.at(6) == middleProcessor2);
            BOOST_CHECK(sorted.at(7) == loopFinalizer2);
            BOOST_CHECK(sorted.at(8) == loopInitiator2);
            BOOST_CHECK(sorted.at(9) == middleProcessor2);
            BOOST_CHECK(sorted.at(10) == loopFinalizer2);
        BOOST_CHECK(sorted.at(11) == middleProcessor);
        BOOST_CHECK(sorted.at(12) == loopFinalizer);
        BOOST_CHECK(sorted.at(13) == loopInitiator);
            // inner loop
            BOOST_CHECK(sorted.at(14) == loopInitiator2);
            BOOST_CHECK(sorted.at(15) == middleProcessor2);
            BOOST_CHECK(sorted.at(16) == loopFinalizer2);
            BOOST_CHECK(sorted.at(17) == loopInitiator2);
            BOOST_CHECK(sorted.at(18) == middleProcessor2);
            BOOST_CHECK(sorted.at(19) == loopFinalizer2);
            BOOST_CHECK(sorted.at(20) == loopInitiator2);
            BOOST_CHECK(sorted.at(21) == middleProcessor2);
            BOOST_CHECK(sorted.at(22) == loopFinalizer2);
        BOOST_CHECK(sorted.at(23) == middleProcessor);
        BOOST_CHECK(sorted.at(24) == loopFinalizer);
    BOOST_CHECK(sorted.at(25) == endProcessor);
}

// Network = (StartProcessor->{LoopInitiator->(LoopInitiator2->MiddleProcessor2->LoopFinalizer2}x3->MiddleProcessor->LoopFinalizer)x2->EndProcessor)
// NOTE: LoopInitiator is connected to LoopFinalizer2, LoopInitiator2 is connected to LoopFinalizer => crossed loops should not be unrolled
BOOST_AUTO_TEST_CASE(crossedLoops)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(loopInitiator);
    processors.push_back(loopInitiator2);
    processors.push_back(middleProcessor);
    processors.push_back(middleProcessor2);
    processors.push_back(startProcessor);
    processors.push_back(loopFinalizer);
    processors.push_back(loopFinalizer2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, loopInitiator2->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator2->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, loopFinalizer2->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer2->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));

    BOOST_CHECK(network->connectPorts(loopFinalizer2->loopOutport_, loopInitiator->loopInport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutport_, loopInitiator2->loopInport_));

    loopInitiator2->loopInport_.setNumLoopIterations(3);
    loopInitiator->loopInport_.setNumLoopIterations(2);

    NetworkGraph graph(processors);

    //BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(loopInitiator));
    BOOST_CHECK(graph.containsProcessor(loopInitiator2));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(middleProcessor2));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer2));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // check predecessors/successors
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(loopInitiator)    &&
                graph.getPredecessors(endProcessor).count(loopInitiator2)   &&
                graph.getPredecessors(endProcessor).count(middleProcessor)  &&
                graph.getPredecessors(endProcessor).count(middleProcessor2) &&
                graph.getPredecessors(endProcessor).count(loopFinalizer)    &&
                graph.getPredecessors(endProcessor).count(loopFinalizer2)   &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );
    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(loopInitiator)    &&
                graph.getSuccessors(startProcessor).count(loopInitiator2)   &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor).count(middleProcessor2) &&
                graph.getSuccessors(startProcessor).count(loopFinalizer)    &&
                graph.getSuccessors(startProcessor).count(loopFinalizer2)   &&
                graph.getSuccessors(startProcessor).count(endProcessor)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    // no loop unrolling = just linear pipeline
    BOOST_CHECK(graph.isSuccessor(loopInitiator, loopFinalizer));
    BOOST_CHECK(graph.isSuccessor(loopInitiator2, loopFinalizer2));
    BOOST_CHECK(graph.isSuccessor(loopInitiator, loopInitiator2));
    BOOST_CHECK(!graph.isSuccessor(loopInitiator2, loopInitiator));
    BOOST_CHECK(!graph.isSuccessor(loopFinalizer, loopInitiator));
    BOOST_CHECK(!graph.isSuccessor(loopFinalizer2, loopInitiator2));
    BOOST_CHECK(!graph.isSuccessor(loopFinalizer2, loopInitiator));
    BOOST_CHECK(graph.isPathElement(middleProcessor, loopInitiator, loopFinalizer));
    BOOST_CHECK(graph.isPathElement(middleProcessor2, loopInitiator2, loopFinalizer2));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, loopInitiator, loopInitiator2));
    BOOST_CHECK(!graph.isPathElement(middleProcessor, loopInitiator2, loopInitiator));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, loopFinalizer, loopInitiator));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, loopFinalizer2, loopInitiator2));
    BOOST_CHECK(!graph.isPathElement(middleProcessor2, loopFinalizer2, loopInitiator));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 8);
    BOOST_CHECK(sorted.at(0) == startProcessor);
    BOOST_CHECK(sorted.at(1) == loopInitiator);
    BOOST_CHECK(sorted.at(2) == loopInitiator2);
    BOOST_CHECK(sorted.at(3) == middleProcessor2);
    BOOST_CHECK(sorted.at(4) == loopFinalizer2);
    BOOST_CHECK(sorted.at(5) == middleProcessor);
    BOOST_CHECK(sorted.at(6) == loopFinalizer);
    BOOST_CHECK(sorted.at(7) == endProcessor);
}

// Network = (StartProcessor->(LoopInitiator->MiddleProcessor->LoopFinalizer)x2->EndProcessor)
//           (StartProcessor2------------------^                                             )
BOOST_AUTO_TEST_CASE(loopExternalDataInflow)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor);
    processors.push_back(loopInitiator);
    processors.push_back(middleProcessor);
    processors.push_back(startProcessor);
    processors.push_back(loopFinalizer);
    processors.push_back(startProcessor2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(startProcessor2->outport_, middleProcessor->inport2_));

    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutport_, loopInitiator->loopInport_));
    loopInitiator->loopInport_.setNumLoopIterations(2);

    NetworkGraph graph(processors);

    //BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(startProcessor2));
    BOOST_CHECK(graph.containsProcessor(loopInitiator));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer));
    BOOST_CHECK(graph.containsProcessor(endProcessor));

    // processors are connected => processors are predecessors/successors of each other
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(startProcessor2)  &&
                graph.getPredecessors(endProcessor).count(loopInitiator)    &&
                graph.getPredecessors(endProcessor).count(middleProcessor)  &&
                graph.getPredecessors(endProcessor).count(loopFinalizer)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                !graph.getSuccessors(startProcessor).count(startProcessor2) &&
                graph.getSuccessors(startProcessor).count(loopInitiator)    &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor).count(loopFinalizer)    &&
                graph.getSuccessors(startProcessor).count(endProcessor)
        );
    BOOST_CHECK(!graph.getSuccessors(startProcessor2).count(startProcessor) &&
                graph.getSuccessors(startProcessor2).count(startProcessor2) &&
                graph.getSuccessors(startProcessor2).count(loopInitiator)   &&
                graph.getSuccessors(startProcessor2).count(middleProcessor) &&
                graph.getSuccessors(startProcessor2).count(loopFinalizer)   &&
                graph.getSuccessors(startProcessor2).count(endProcessor)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor2, middleProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor2, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor2, loopInitiator));                     //< loop unrolling
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor2, loopInitiator));  //< loop unrolling
    BOOST_CHECK(!graph.isPathElement(middleProcessor, loopInitiator, startProcessor2));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 9);
    BOOST_CHECK(beforeInVector(sorted, startProcessor2, middleProcessor));

    sorted.erase(std::find(sorted.begin(), sorted.end(), startProcessor2)); //< remove startprocessor2 before checking loop
    BOOST_CHECK(sorted.at(0) == startProcessor);
    BOOST_CHECK(sorted.at(1) == loopInitiator);
    BOOST_CHECK(sorted.at(2) == middleProcessor);
    BOOST_CHECK(sorted.at(3) == loopFinalizer);
    BOOST_CHECK(sorted.at(4) == loopInitiator);
    BOOST_CHECK(sorted.at(5) == middleProcessor);
    BOOST_CHECK(sorted.at(6) == loopFinalizer);
    BOOST_CHECK(sorted.at(7) == endProcessor);
}

// Network = (StartProcessor->(LoopInitiator->MiddleProcessor->LoopFinalizer)x2->EndProcessor )
//           (                                             \-------------------->EndProcessor2)
BOOST_AUTO_TEST_CASE(loopExternalDataOutflow)
{
    vector<Processor*> processors;
    processors.push_back(endProcessor2);
    processors.push_back(endProcessor);
    processors.push_back(loopInitiator);
    processors.push_back(middleProcessor);
    processors.push_back(startProcessor);
    processors.push_back(loopFinalizer);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport2_, endProcessor2->inport_));

    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutport_, loopInitiator->loopInport_));
    loopInitiator->loopInport_.setNumLoopIterations(2);

    NetworkGraph graph(processors);

    //BOOST_CHECK(graph.getProcessors().size() == 4);
    BOOST_CHECK(graph.containsProcessor(startProcessor));
    BOOST_CHECK(graph.containsProcessor(loopInitiator));
    BOOST_CHECK(graph.containsProcessor(middleProcessor));
    BOOST_CHECK(graph.containsProcessor(loopFinalizer));
    BOOST_CHECK(graph.containsProcessor(endProcessor));
    BOOST_CHECK(graph.containsProcessor(endProcessor2));

    // check predecessor/successors
    BOOST_CHECK(graph.getPredecessors(endProcessor).count(startProcessor)   &&
                graph.getPredecessors(endProcessor).count(loopInitiator)    &&
                graph.getPredecessors(endProcessor).count(middleProcessor)  &&
                graph.getPredecessors(endProcessor).count(loopFinalizer)    &&
                graph.getPredecessors(endProcessor).count(endProcessor)     &&
                !graph.getPredecessors(endProcessor).count(endProcessor2)
        );
    BOOST_CHECK(graph.getPredecessors(endProcessor2).count(startProcessor)  &&
                graph.getPredecessors(endProcessor2).count(loopInitiator)   &&
                graph.getPredecessors(endProcessor2).count(middleProcessor) &&
                graph.getPredecessors(endProcessor2).count(loopFinalizer)   &&
                !graph.getPredecessors(endProcessor2).count(endProcessor)   &&
                graph.getPredecessors(endProcessor2).count(endProcessor2)
        );

    BOOST_CHECK(graph.getSuccessors(startProcessor).count(startProcessor)   &&
                graph.getSuccessors(startProcessor).count(loopInitiator)    &&
                graph.getSuccessors(startProcessor).count(middleProcessor)  &&
                graph.getSuccessors(startProcessor).count(loopFinalizer)    &&
                graph.getSuccessors(startProcessor).count(endProcessor)     &&
                graph.getSuccessors(startProcessor).count(endProcessor2)
        );

    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor));
    BOOST_CHECK(graph.isSuccessor(startProcessor, endProcessor2));
    BOOST_CHECK(graph.isSuccessor(startProcessor, middleProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor));
    BOOST_CHECK(graph.isPathElement(middleProcessor, startProcessor, endProcessor2));
    BOOST_CHECK(graph.isSuccessor(middleProcessor, loopInitiator));                     //< loop unrolling
    BOOST_CHECK(graph.isPathElement(loopFinalizer, middleProcessor, loopInitiator));    //< loop unrolling
    BOOST_CHECK(!graph.isPathElement(endProcessor, middleProcessor, loopInitiator));

    // check topological sorting
    std::vector<Processor*> sorted = graph.sortTopologically();
    BOOST_CHECK(sorted.size() == 9);
    BOOST_CHECK(beforeInVector(sorted, middleProcessor, endProcessor2));

    sorted.erase(std::find(sorted.begin(), sorted.end(), endProcessor2)); //< remove endprocessor2 before checking loop
    BOOST_CHECK(sorted.at(0) == startProcessor);
    BOOST_CHECK(sorted.at(1) == loopInitiator);
    BOOST_CHECK(sorted.at(2) == middleProcessor);
    BOOST_CHECK(sorted.at(3) == loopFinalizer);
    BOOST_CHECK(sorted.at(4) == loopInitiator);
    BOOST_CHECK(sorted.at(5) == middleProcessor);
    BOOST_CHECK(sorted.at(6) == loopFinalizer);
    BOOST_CHECK(sorted.at(7) == endProcessor);
}

BOOST_AUTO_TEST_SUITE_END()


//-------------------------------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(NetworkEvaluatorTests, TestFixture)

// Checks that processors are initialized by the evaluator
BOOST_AUTO_TEST_CASE(networkInitialize)
{
    BOOST_CHECK(!startProcessor->isInitialized());
    BOOST_CHECK(!startProcessor2->isInitialized());
    BOOST_CHECK(!middleProcessor->isInitialized());
    BOOST_CHECK(!middleProcessor2->isInitialized());
    BOOST_CHECK(!endProcessor->isInitialized());
    BOOST_CHECK(!endProcessor2->isInitialized());
    BOOST_CHECK(!loopInitiator->isInitialized());
    BOOST_CHECK(!loopInitiator2->isInitialized());
    BOOST_CHECK(!loopFinalizer->isInitialized());
    BOOST_CHECK(!loopFinalizer2->isInitialized());

    evaluator->setProcessorNetwork(network);

    BOOST_CHECK(startProcessor->isInitialized());
    BOOST_CHECK(startProcessor2->isInitialized());
    BOOST_CHECK(middleProcessor->isInitialized());
    BOOST_CHECK(middleProcessor2->isInitialized());
    BOOST_CHECK(endProcessor->isInitialized());
    BOOST_CHECK(endProcessor2->isInitialized());
    BOOST_CHECK(loopInitiator->isInitialized());
    BOOST_CHECK(loopInitiator2->isInitialized());
    BOOST_CHECK(loopFinalizer->isInitialized());
    BOOST_CHECK(loopFinalizer2->isInitialized());
}

// Checks that processors are deinitialized by the evaluator on network removal
BOOST_AUTO_TEST_CASE(networkDeinitialize)
{
    evaluator->setProcessorNetwork(network);
    evaluator->setProcessorNetwork(0, true);

    BOOST_CHECK(!startProcessor->isInitialized());
    BOOST_CHECK(!startProcessor2->isInitialized());
    BOOST_CHECK(!middleProcessor->isInitialized());
    BOOST_CHECK(!middleProcessor2->isInitialized());
    BOOST_CHECK(!endProcessor->isInitialized());
    BOOST_CHECK(!endProcessor2->isInitialized());
    BOOST_CHECK(!loopInitiator->isInitialized());
    BOOST_CHECK(!loopInitiator2->isInitialized());
    BOOST_CHECK(!loopFinalizer->isInitialized());
    BOOST_CHECK(!loopFinalizer2->isInitialized());
}

// Checks evaluation order of an unconnected network
BOOST_AUTO_TEST_CASE(unconnectedNetworkEvaluation)
{
    evaluator->setProcessorNetwork(network);
    BOOST_CHECK(evalOrderRec->evalOrder.empty());

    evaluator->process();

    // expected:
    // - only end processors have been evaluated
    BOOST_CHECK(evalOrderRec->evalOrder.size() == 2);
    BOOST_CHECK(inVector(evalOrderRec->evalOrder, endProcessor));
    BOOST_CHECK(inVector(evalOrderRec->evalOrder, endProcessor2));
}

// Checks evaluation order of a linear pipeline
// Network = (StartProcessor->MiddleProcessor->MiddleProcessor2->EndProcessor;EndProcessor2)
BOOST_AUTO_TEST_CASE(linearPipelineEvaluation)
{
    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, endProcessor->inport_));

    evaluator->setProcessorNetwork(network);
    BOOST_CHECK(evalOrderRec->evalOrder.empty());

    evaluator->process();

    // expected network evaluation:
    // - only processors in the pipeline and the endProcessor2 have been evaluated
    // - evaluated processors are valid
    BOOST_CHECK(evalOrderRec->evalOrder.size() == 5);
    BOOST_CHECK(beforeInVector(evalOrderRec->evalOrder, startProcessor, middleProcessor));
    BOOST_CHECK(beforeInVector(evalOrderRec->evalOrder, middleProcessor, middleProcessor2));
    BOOST_CHECK(beforeInVector(evalOrderRec->evalOrder, middleProcessor2, endProcessor));
    BOOST_CHECK(inVector(evalOrderRec->evalOrder, endProcessor2));

    BOOST_CHECK(startProcessor->isValid());
    BOOST_CHECK(middleProcessor->isValid());
    BOOST_CHECK(middleProcessor2->isValid());
    BOOST_CHECK(endProcessor->isValid());
    BOOST_CHECK(endProcessor2->isValid());
}

// Checks evaluation order of a simple (linear, non-nested) loop
// Network = (StartProcessor->(LoopInitiator->MiddleProcessor->LoopFinalizer)x3->EndProcessor)
BOOST_AUTO_TEST_CASE(simpleLoopEvaluation)
{
    // remove endProcessor2 from network for being able to precisely check the evaluation order
    network->removeProcessor(endProcessor2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutport_, loopInitiator->loopInport_));
    loopInitiator->loopInport_.setNumLoopIterations(3);

    evaluator->setProcessorNetwork(network);
    BOOST_CHECK(evalOrderRec->evalOrder.empty());

    evaluator->process();

    // check evaluation order
    BOOST_CHECK(evalOrderRec->evalOrder.size() == 11);
    BOOST_CHECK(evalOrderRec->loopIterations.size() == 11);
    BOOST_CHECK(evalOrderRec->evalOrder.at(0) == startProcessor);

    BOOST_CHECK(evalOrderRec->evalOrder.at(1) == loopInitiator);
    BOOST_CHECK(evalOrderRec->loopIterations.at(1) == 0);
    BOOST_CHECK(evalOrderRec->evalOrder.at(2) == middleProcessor);
    BOOST_CHECK(evalOrderRec->evalOrder.at(3) == loopFinalizer);
    BOOST_CHECK(evalOrderRec->loopIterations.at(3) == 0);

    BOOST_CHECK(evalOrderRec->evalOrder.at(4) == loopInitiator);
    BOOST_CHECK(evalOrderRec->loopIterations.at(4) == 1);
    BOOST_CHECK(evalOrderRec->evalOrder.at(5) == middleProcessor);
    BOOST_CHECK(evalOrderRec->evalOrder.at(6) == loopFinalizer);
    BOOST_CHECK(evalOrderRec->loopIterations.at(6) == 1);

    BOOST_CHECK(evalOrderRec->evalOrder.at(7) == loopInitiator);
    BOOST_CHECK(evalOrderRec->loopIterations.at(7) == 2);
    BOOST_CHECK(evalOrderRec->evalOrder.at(8) == middleProcessor);
    BOOST_CHECK(evalOrderRec->evalOrder.at(9) == loopFinalizer);
    BOOST_CHECK(evalOrderRec->loopIterations.at(9) == 2);

    BOOST_CHECK(evalOrderRec->evalOrder.at(10) == endProcessor);

    // evaluated processors are expected to be valid, except loopInitiator (invalidated by loopFinalizer)
    BOOST_CHECK(startProcessor->isValid());
    BOOST_CHECK(!loopInitiator->isValid());
    BOOST_CHECK(middleProcessor->isValid());
    BOOST_CHECK(loopFinalizer->isValid());
    BOOST_CHECK(endProcessor->isValid());
}

// Checks evaluation order of a nested loop
// Network = (StartProcessor->(LoopInitiator->(LoopInitiator2->MiddleProcessor2->LoopFinalizer2)x3->MiddleProcessor->LoopFinalizer)x2->EndProcessor)
BOOST_AUTO_TEST_CASE(nestedLoopEvaluation)
{
    // remove endProcessor2 from network for being able to precisely check the evaluation order
    network->removeProcessor(endProcessor2);

    BOOST_CHECK(network->connectPorts(startProcessor->outport_, loopInitiator->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator->outport_, loopInitiator2->inport_));
    BOOST_CHECK(network->connectPorts(loopInitiator2->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, loopFinalizer2->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer2->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, loopFinalizer->inport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->outport_, endProcessor->inport_));

    BOOST_CHECK(network->connectPorts(loopFinalizer2->loopOutport_, loopInitiator2->loopInport_));
    BOOST_CHECK(network->connectPorts(loopFinalizer->loopOutport_, loopInitiator->loopInport_));

    loopInitiator2->loopInport_.setNumLoopIterations(3);
    loopInitiator->loopInport_.setNumLoopIterations(2);

    evaluator->setProcessorNetwork(network);
    BOOST_CHECK(evalOrderRec->evalOrder.empty());

    evaluator->process();

    // check evaluation order
    BOOST_CHECK(evalOrderRec->evalOrder.size() == 26);
    BOOST_CHECK(evalOrderRec->loopIterations.size() == 26);
    BOOST_CHECK(evalOrderRec->evalOrder.at(0) == startProcessor);
    // outer loop
    BOOST_CHECK(evalOrderRec->evalOrder.at(1) == loopInitiator);
    BOOST_CHECK(evalOrderRec->loopIterations.at(1) == 0);

        // inner loop
        BOOST_CHECK(evalOrderRec->evalOrder.at(2) == loopInitiator2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(2) == 0);
        BOOST_CHECK(evalOrderRec->evalOrder.at(3) == middleProcessor2);
        BOOST_CHECK(evalOrderRec->evalOrder.at(4) == loopFinalizer2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(4) == 0);
        BOOST_CHECK(evalOrderRec->evalOrder.at(5) == loopInitiator2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(5) == 1);
        BOOST_CHECK(evalOrderRec->evalOrder.at(6) == middleProcessor2);
        BOOST_CHECK(evalOrderRec->evalOrder.at(7) == loopFinalizer2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(7) == 1);
        BOOST_CHECK(evalOrderRec->evalOrder.at(8) == loopInitiator2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(8) == 2);
        BOOST_CHECK(evalOrderRec->evalOrder.at(9) == middleProcessor2);
        BOOST_CHECK(evalOrderRec->evalOrder.at(10) == loopFinalizer2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(10) == 2);

    BOOST_CHECK(evalOrderRec->evalOrder.at(11) == middleProcessor);
    BOOST_CHECK(evalOrderRec->evalOrder.at(12) == loopFinalizer);
    BOOST_CHECK(evalOrderRec->loopIterations.at(12) == 0);
    BOOST_CHECK(evalOrderRec->evalOrder.at(13) == loopInitiator);
    BOOST_CHECK(evalOrderRec->loopIterations.at(13) == 1);

        // inner loop
        BOOST_CHECK(evalOrderRec->evalOrder.at(14) == loopInitiator2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(14) == 0);
        BOOST_CHECK(evalOrderRec->evalOrder.at(15) == middleProcessor2);
        BOOST_CHECK(evalOrderRec->evalOrder.at(16) == loopFinalizer2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(16) == 0);
        BOOST_CHECK(evalOrderRec->evalOrder.at(17) == loopInitiator2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(17) == 1);
        BOOST_CHECK(evalOrderRec->evalOrder.at(18) == middleProcessor2);
        BOOST_CHECK(evalOrderRec->evalOrder.at(19) == loopFinalizer2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(19) == 1);
        BOOST_CHECK(evalOrderRec->evalOrder.at(20) == loopInitiator2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(20) == 2);
        BOOST_CHECK(evalOrderRec->evalOrder.at(21) == middleProcessor2);
        BOOST_CHECK(evalOrderRec->evalOrder.at(22) == loopFinalizer2);
        BOOST_CHECK(evalOrderRec->loopIterations.at(22) == 2);

    BOOST_CHECK(evalOrderRec->evalOrder.at(23) == middleProcessor);
    BOOST_CHECK(evalOrderRec->evalOrder.at(24) == loopFinalizer);
    BOOST_CHECK(evalOrderRec->loopIterations.at(24) == 1);

    BOOST_CHECK(evalOrderRec->evalOrder.at(25) == endProcessor);

    // evaluated processors are expected to be valid, except loopInitiators (invalidated by loopFinalizers)
    BOOST_CHECK(startProcessor->isValid());
    BOOST_CHECK(!loopInitiator->isValid());
    BOOST_CHECK(!loopInitiator2->isValid());
    BOOST_CHECK(middleProcessor2->isValid());
    BOOST_CHECK(loopFinalizer2->isValid());
    BOOST_CHECK(middleProcessor->isValid());
    BOOST_CHECK(loopFinalizer->isValid());
    BOOST_CHECK(endProcessor->isValid());
}

// Checks evaluation of linear pipeline with valid processors (which should not be processed)
// Network = (StartProcessor->MiddleProcessor->MiddleProcessor2->EndProcessor;EndProcessor2)
BOOST_AUTO_TEST_CASE(validProcessorEvaluation)
{
    BOOST_CHECK(network->connectPorts(startProcessor->outport_, middleProcessor->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor->outport_, middleProcessor2->inport_));
    BOOST_CHECK(network->connectPorts(middleProcessor2->outport_, endProcessor->inport_));

    evaluator->setProcessorNetwork(network);
    BOOST_CHECK(evalOrderRec->evalOrder.empty());

    // set start and middle processor valid
    startProcessor->setValid();
    middleProcessor->setValid();

    evaluator->process();

    // expected network evaluation:
    // - only invalid! processors in the pipeline and the endProcessor2 have been evaluated
    // - evaluated processors are valid
    BOOST_CHECK(evalOrderRec->evalOrder.size() == 3);
    BOOST_CHECK(beforeInVector(evalOrderRec->evalOrder, middleProcessor2, endProcessor));
    BOOST_CHECK(inVector(evalOrderRec->evalOrder, endProcessor2));

    BOOST_CHECK(startProcessor->isValid());
    BOOST_CHECK(middleProcessor->isValid());
    BOOST_CHECK(middleProcessor2->isValid());
    BOOST_CHECK(endProcessor->isValid());
    BOOST_CHECK(endProcessor2->isValid());
}

BOOST_AUTO_TEST_SUITE_END()

