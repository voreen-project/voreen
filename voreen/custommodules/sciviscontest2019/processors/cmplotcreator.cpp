#include "cmplotcreator.h"

#include <algorithm>

namespace voreen {


	CMPlotCreator::CMPlotCreator()
		: AsyncComputeProcessor()
		, inport_(Port::INPORT, "cmplotdata", "Particle Data Input")
		, outport_(Port::OUTPORT, "outport", "Particle Plot Output")
	{
		addPort(inport_);
		addPort(outport_);


		std::vector<float> sampleData = std::vector<float>();
		float max = 0.0f;
		for (int i = 0; i < 150; i++) {
			float dataPoint = rand() % 2000;
			sampleData.push_back(dataPoint);
			max = std::max(dataPoint, max);
		}
		data = { "Data1", sampleData, 0.0f, max };
	}

	CMPlotCreator::~CMPlotCreator()
	{
	}

	Processor* CMPlotCreator::create() const
	{
		return new CMPlotCreator();
	}

	CMParticleData* CMPlotCreator::prepareComputeInput()
	{
		//const CMParticleData* inputPtr = inport_.getThreadSafeData();
		const CMParticleData* inputPtr = inport_.getData();
		if (!inputPtr)
			throw InvalidInputException("No input", InvalidInputException::S_WARNING);

		return (CMParticleData*)inputPtr;
	}

	CMPlotData* CMPlotCreator::compute(CMParticleData* input, ProgressReporter & progressReporter) const
	{
		CMPlotData* dataOut = new CMPlotData();

		std::vector<float> starVector = std::vector<float>();
		std::vector<float> windVector = std::vector<float>();
		std::vector<float> gasVector = std::vector<float>();
		std::vector<float> agnVector = std::vector<float>();

		std::vector<float> multiVector = std::vector<float>();

		std::vector<float> starNewVector = std::vector<float>();
		std::vector<float> windNewVector = std::vector<float>();
		std::vector<float> gasNewVector = std::vector<float>();
		std::vector<float> agnNewVector = std::vector<float>();

		std::vector<int64_t> starIDs = std::vector<int64_t>();
		std::vector<int64_t> windIDs = std::vector<int64_t>();
		std::vector<int64_t> gasIDs = std::vector<int64_t>();
		std::vector<int64_t> agnIDs = std::vector<int64_t>();

		const char* enabledState = input->getEnabledState();

		for (int i = 0; i < 625; i++) {
			CMParticleDataTimeSlice* slice = input->sliceAtTimeStep(i);
			const CMParticle* particleArray = slice->startUsingParticles();

			int numStar = 0;
			int numWind = 0;
			int numGas = 0;
			int numAgn = 0;

			int numMulti = 0;
			
			std::vector<int64_t> newStarIDs = std::vector<int64_t>();
			std::vector<int64_t> newWindIDs = std::vector<int64_t>();
			std::vector<int64_t> newGasIDs = std::vector<int64_t>();
			std::vector<int64_t> newAgnIDs = std::vector<int64_t>();

			for (int x = 0; x < slice->getNumberOfParticles(); x++) {
				if (!enabledState[particleArray[x].ident])
					continue;

				uint16_t mask = particleArray[x].mask;
				if (mask & 2) {
					int multi = 0;
					if (mask & 32) {
						numStar++;
						multi++;
						newStarIDs.push_back(particleArray[x].ident);
					}
					if (mask & 64) {
						numWind++;
						multi++;
						newWindIDs.push_back(particleArray[x].ident);
					}
					if (mask & 128) {
						numGas++;
						multi++;
						newGasIDs.push_back(particleArray[x].ident);
					}
					if (multi > 1) {
						numMulti++;
					}
				}
				else {
					if (mask & 256) {
						numAgn++;
						newAgnIDs.push_back(particleArray[x].ident);
					}
				}

			}

			std::vector<int64_t> starDiff = std::vector<int64_t>(std::max(starIDs.size(), newStarIDs.size()));
			std::sort(newStarIDs.begin(), newStarIDs.end());
			std::set_difference(newStarIDs.begin(), newStarIDs.end(), starIDs.begin(), starIDs.end(), starDiff.begin());

			std::vector<int64_t> windDiff = std::vector<int64_t>(std::max(windIDs.size(), newWindIDs.size()));
			std::sort(newWindIDs.begin(), newWindIDs.end());
			std::set_difference(newWindIDs.begin(), newWindIDs.end(), windIDs.begin(), windIDs.end(), windDiff.begin());

			std::vector<int64_t> gasDiff = std::vector<int64_t>(std::max(gasIDs.size(), newGasIDs.size()));
			std::sort(newGasIDs.begin(), newGasIDs.end());
			std::set_difference(newGasIDs.begin(), newGasIDs.end(), gasIDs.begin(), gasIDs.end(), gasDiff.begin());

			std::vector<int64_t> agnDiff = std::vector<int64_t>(std::max(agnIDs.size(), newAgnIDs.size()));
			std::sort(newAgnIDs.begin(), newAgnIDs.end());
			std::set_difference(newAgnIDs.begin(), newAgnIDs.end(), agnIDs.begin(), agnIDs.end(), agnDiff.begin());

			starIDs = newStarIDs;
			windIDs = newWindIDs;
			gasIDs = newGasIDs;
			agnIDs = newAgnIDs;

			starNewVector.push_back((float)count_if(starDiff.begin(), starDiff.end(), [](int64_t i) { return i != 0; }));
			windNewVector.push_back((float)count_if(windDiff.begin(), windDiff.end(), [](int64_t i) { return i != 0; }));
			gasNewVector.push_back((float)count_if(gasDiff.begin(), gasDiff.end(), [](int64_t i) { return i != 0; }));
			agnNewVector.push_back((float)count_if(agnDiff.begin(), agnDiff.end(), [](int64_t i) { return i != 0; }));

			starVector.push_back((float)numStar);
			windVector.push_back((float)numWind);
			gasVector.push_back((float)numGas);
			agnVector.push_back((float)numAgn);

			multiVector.push_back((float)numMulti);

			slice->finishedUsingParticles();

			progressReporter.setProgress(i / 625.0f);
		}

		float maxStar = 0.0f;
		float maxWind = 0.0f;
		float maxGas = 0.0f;
		float maxAgn = 0.0f;

		float maxMulti = 0.0f;

		float maxNewStar = 0.0f;
		float maxNewWind = 0.0f;
		float maxNewGas = 0.0f;
		float maxNewAgn = 0.0f;

		for (int i = 0; i < starVector.size(); i++) {
			maxStar = std::max(maxStar, starVector[i]);
			maxWind = std::max(maxWind, windVector[i]);
			maxGas = std::max(maxGas, gasVector[i]);
			maxAgn = std::max(maxAgn, agnVector[i]);
			
			maxMulti = std::max(maxMulti, multiVector[i]);

			maxNewStar = std::max(maxNewStar, starNewVector[i]);
			maxNewWind = std::max(maxNewWind, windNewVector[i]);
			maxNewGas = std::max(maxNewGas, gasNewVector[i]);
			maxNewAgn = std::max(maxNewAgn, agnNewVector[i]);
		}

		CMPlotDataRow starRow = { "Number of Star Particles", starVector, 0.0f, maxStar };
		CMPlotDataRow windRow = { "Number of Wind Particles", windVector, 0.0f, maxWind };
		CMPlotDataRow gasRow = { "Number of Gas Particles", gasVector, 0.0f, maxGas };
		CMPlotDataRow agnRow = { "Number of AGN Dark Mater", agnVector, 0.0f, maxAgn };

		CMPlotDataRow multiRow = { "Number of Particles in Multiple Categories", multiVector, 0.0f, maxMulti };

		CMPlotDataRow starNewRow = { "Number of new Star Particles", starNewVector, 0.0f, maxNewStar };
		CMPlotDataRow windNewRow = { "Number of new Wind Particles", windNewVector, 0.0f, maxNewWind };
		CMPlotDataRow gasNewRow = { "Number of new Gas Particles", gasNewVector, 0.0f, maxNewGas };
		CMPlotDataRow agnNewRow = { "Number of new AGN Dark Mater", agnNewVector, 0.0f, maxNewAgn };

		dataOut->addDataRow(starRow);
		dataOut->addDataRow(windRow);
		dataOut->addDataRow(gasRow);
		dataOut->addDataRow(agnRow);

		dataOut->addDataRow(multiRow);

		dataOut->addDataRow(starNewRow);
		dataOut->addDataRow(windNewRow);
		dataOut->addDataRow(gasNewRow);
		dataOut->addDataRow(agnNewRow);

		return dataOut;
	}

	void CMPlotCreator::processComputeOutput(CMPlotData* output)
	{
		outport_.setData(output, true);
	}

}