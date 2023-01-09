import voreen

ensemble_path = '/DATA/groupdrive/scivis_contest/2022/H5Reparametrized/'
output_name = 'similarity.vsm'

voreen.setPropertyValue('EnsembleDataSource', 'ensemblepath', ensemble_path)
voreen.setPropertyValue('EnsembleDataSource', 'loadDataset', 1)
voreen.render()
voreen.setPropertyValue('SimilarityMatrixCreator', 'numSeedPoints', 16384)
voreen.setPropertyValue('SimilarityMatrixCreator', 'manualUpdateButton_', 1)
voreen.render()
voreen.setPropertyValue('SimilarityMatrixSave', 'filenameprop', ensemble_path + output_name)
voreen.setPropertyValue('SimilarityMatrixSave', 'saveButton', 1)
voreen.render()
