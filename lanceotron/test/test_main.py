import unittest
import lanceotron.utils
import pkg_resources
import os
from shutil import rmtree

import numpy as np
import pandas as pd

from lanceotron import find_and_score_peaks
from tensorflow import keras

# Simple test run
simple_run = {
    "file": "test/chr22.bw",
    "folder": "./",
    "threshold": 4,
    "window": 400,
    "skipheader": False,
}

no_dir_run = {
    "file": "test/chr22.bw",
    "folder": "./test/not_a_folder",
    "threshold": 4,
    "window": 400,
    "skipheader": False, 
}

# Simple run through on toy data
class TestAll(unittest.TestCase):
    """Full run through on toy data"""

    def test_example(self):
        find_and_score_peaks(**simple_run)

    def test_new_folder(self):
        find_and_score_peaks(**no_dir_run)
        rmtree(no_dir_run['folder'])


class TestOutput(unittest.TestCase):
    # Little hack for avoiding running for each test case
    @classmethod
    def setUpClass(cls):
        super(TestOutput, cls).setUpClass()
        find_and_score_peaks(**simple_run)
        cls.output = pd.read_csv(
            simple_run["file"].split("/")[-1].replace(".bw", "_L-tron.bed"), sep="\t"
        )

    def test_colnames(self):
        for col in ["chrom", "start", "end", "overall_peak_score", "shape_score"]:
            assert col in self.output.columns

    def test_nrow(self):
        assert self.output.shape[0] > 1

    def test_ncol(self):
        assert self.output.shape[1] == 17

    def test_coltypes(self):
        assert str(self.output["chrom"].dtype) == "object"
        assert str(self.output["start"].dtype) == "int64"
        assert str(self.output["end"].dtype) == "int64"
        assert str(self.output["shape_score"].dtype) == "float64"
    
    def test_output(self):
        assert os.path.isdir(lanceotron.utils.make_directory_name('test/test_directory'))
        rmtree('test/test_directory')


class TestResources(unittest.TestCase):
    """Ensure that all of the resources are available."""

    def test_model(self):
        assert os.path.isfile(
            pkg_resources.resource_filename(
                "lanceotron.static", "wide_and_deep_fully_trained_v5_03.h5"
            )
        )

    def test_widescaler(self):
        assert os.path.isfile(
            pkg_resources.resource_filename(
                "lanceotron.static", "standard_scaler_wide_v5_03.p"
            )
        )

    def test_deepscaler(self):
        assert os.path.isfile(
            pkg_resources.resource_filename(
                "lanceotron.static", "standard_scaler_deep_v5_03.p"
            )
        )


class TestModel(unittest.TestCase):
    """Ensure that the model is correctly formed"""

    def setUp(self):
        self.model = lanceotron.utils.build_model()

    def test_nlayer(self):
        assert (
            len(self.model.layers) == 37
        ), "Did the model change? It has the wrong number of layers."

    def test_lastlayer(self):
        assert isinstance(
            self.model.layers[-1], keras.layers.Dense
        ), "The last layer is not a dense prediction layer. Did the model load correctly?"

    def test_hasweights(self):
        np.testing.assert_almost_equal(
            self.model.layers[-1].get_weights()[1][1],
            -0.015733806,
            err_msg="The model does not have the correct weights. Check the weight file.",
        )