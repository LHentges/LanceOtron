import numpy as np
from scipy.stats import poisson
import pkg_resources
from pathlib import Path


def make_directory_name(directory: str) -> str:
    """Ensures a trailing slash on directory names.

    Args:
        directory (str): Path to the directory

    Returns:
        str: Directory name with a guarenteed trailing slash
    """
    if directory[-1] != "/":
        directory += "/"

    # Create the directory if it doens't exist
    Path(directory).mkdir(parents = True, exist_ok = True)

    return directory


def build_model():
    """
    Build and return the Lanceotron Keras model.
    """
    from tensorflow import keras
    import tensorflow.keras.backend as K

    deep_dense_size = 10
    dropout_rate = 0.5
    first_filter_num = 70
    first_filter_size = 9
    hidden_filter_num = 120
    hidden_filter_size = 6
    learning_rate = 0.0001
    wide_and_deep_dense_size = 70

    input_wide = keras.layers.Input(shape=(12, 1))
    wide_model = keras.layers.Flatten()(input_wide)
    input_deep = keras.layers.Input((2000, 1))
    deep_model = input_deep

    # deep model first conv layer
    deep_model = keras.layers.Convolution1D(
        first_filter_num, kernel_size=first_filter_size, padding="same"
    )(deep_model)
    deep_model = keras.layers.BatchNormalization()(deep_model)
    deep_model = keras.layers.LeakyReLU()(deep_model)

    # deep model - 4 conv blocks
    deep_model = keras.layers.Convolution1D(
        hidden_filter_num, kernel_size=hidden_filter_size, padding="same"
    )(deep_model)
    deep_model = keras.layers.BatchNormalization()(deep_model)
    deep_model = keras.layers.LeakyReLU()(deep_model)
    deep_model = keras.layers.MaxPool1D(pool_size=2)(deep_model)

    deep_model = keras.layers.Convolution1D(
        hidden_filter_num, kernel_size=hidden_filter_size, padding="same"
    )(deep_model)
    deep_model = keras.layers.BatchNormalization()(deep_model)
    deep_model = keras.layers.LeakyReLU()(deep_model)
    deep_model = keras.layers.MaxPool1D(pool_size=2)(deep_model)

    deep_model = keras.layers.Convolution1D(
        hidden_filter_num, kernel_size=hidden_filter_size, padding="same"
    )(deep_model)
    deep_model = keras.layers.BatchNormalization()(deep_model)
    deep_model = keras.layers.LeakyReLU()(deep_model)
    deep_model = keras.layers.MaxPool1D(pool_size=2)(deep_model)

    deep_model = keras.layers.Convolution1D(
        hidden_filter_num, kernel_size=hidden_filter_size, padding="same"
    )(deep_model)
    deep_model = keras.layers.BatchNormalization()(deep_model)
    deep_model = keras.layers.LeakyReLU()(deep_model)
    deep_model = keras.layers.MaxPool1D(pool_size=2)(deep_model)

    # deep model - dense layer with dropout
    deep_model = keras.layers.Dense(deep_dense_size)(deep_model)
    deep_model = keras.layers.BatchNormalization()(deep_model)
    deep_model = keras.layers.LeakyReLU()(deep_model)
    deep_model = keras.layers.Dropout(dropout_rate)(deep_model)
    deep_model = keras.layers.Flatten()(deep_model)

    # shape output only dense layer
    shape_output = keras.layers.Dense(
        2, activation="softmax", name="shape_classification"
    )(deep_model)

    # p-value output only dense layer
    pvalue_output = keras.layers.Dense(
        2, activation="softmax", name="pvalue_classification"
    )(wide_model)

    # combine wide and deep paths
    concat = keras.layers.concatenate([wide_model, deep_model, pvalue_output])
    wide_and_deep = keras.layers.Dense(wide_and_deep_dense_size)(concat)
    wide_and_deep = keras.layers.BatchNormalization()(wide_and_deep)
    wide_and_deep = keras.layers.LeakyReLU()(wide_and_deep)
    wide_and_deep = keras.layers.Dense(wide_and_deep_dense_size)(wide_and_deep)
    wide_and_deep = keras.layers.BatchNormalization()(wide_and_deep)
    wide_and_deep = keras.layers.LeakyReLU()(wide_and_deep)
    output = keras.layers.Dense(2, activation="softmax", name="overall_classification")(
        wide_and_deep
    )
    model = keras.models.Model(
        inputs=[input_deep, input_wide], outputs=[output, shape_output, pvalue_output]
    )

    # load model weights
    model.load_weights(
        pkg_resources.resource_filename(
            "lanceotron.static", "wide_and_deep_fully_trained_v5_03.h5"
        )
    )
    return model


def calculate_pvalue_from_input(
    chrom, start, end, seq_depth_test, seq_depth_control, pyBigWig_object, max_height
):
    ave_coverage_input = pyBigWig_object.stats(
        chrom, start, end, type="mean", exact=True
    )[0] * (seq_depth_test / seq_depth_control)
    with np.errstate(divide="ignore"):
        pvalue = -1 * np.log10(1 - poisson.cdf(max_height, ave_coverage_input))
    if np.isinf(pvalue):
        pvalue = 100.0
    if pvalue > 100:
        pvalue = 100.0
    return pvalue
