import numpy as np
from scipy.signal import butter, filtfilt

def readXLSX(filename):
    """
    A function that reads an Excel file and returns the data as a pandas DataFrame.

    Parameters:
        filename (str): The path to the Excel file.

    Returns:
        pandas.DataFrame: The DataFrame containing the data from the Excel file.
    """
    #file_path = os.path.join(os.path.dirname(__file__), 'source', filename)
    df = pd.read_excel(filename, header=0)
    return df



def rescale_linear(array, new_min, new_max, minimum, maximum):
    """
    Rescale an array linearly.

    Args:
        array (numpy.ndarray): The array to be rescaled.
        new_min (float): The new minimum value of the rescaled array.
        new_max (float): The new maximum value of the rescaled array.
        minimum (float): The minimum value of the original array.
        maximum (float): The maximum value of the original array.

    Returns:
        numpy.ndarray: The rescaled array.

    Description:
        This function rescales an array linearly from the range [minimum, maximum] to the range [new_min, new_max].
        The rescaling formula is:
            rescaled_value = m * original_value + b
        where m is the slope and b is the y-intercept.

        The slope is calculated as:
            m = (new_max - new_min) / (maximum - minimum)
        The y-intercept is calculated as:
            b = new_min - m * minimum

        Note: The original array is not modified. A new array is returned with the rescaled values.

    Example:
        >>> array = np.array([1, 2, 3, 4, 5])
        >>> rescaled_array = rescale_linear(array, 0, 10, 1, 5)
        >>> print(rescaled_array)
        [ 0.  2.  4.  6. 10.]
    """
    
    m = (new_max - new_min) / (maximum - minimum)
    b = new_min - m * minimum
    return m * array + b


def smooth_Gaussian(list, degree=5):
    window = degree*2-1
    weight = np.array([1.0]*window)
    weightGauss = []
    for i in range(window):
        frac = (i-degree+1)/float(window)
        gauss = 1/(np.exp((4*(frac))**2))
        weightGauss.append(gauss)
    weight = np.array(weightGauss)*weight
    smoothed = [0.0]*(len(list)-window)
    for i in range(len(smoothed)):
        smoothed[i] = sum(np.array(list[i:i+window])*weight)/sum(weight)

    start = round(window/2)
    end = window - start
    add_start = np.empty(start)
    add_end = np.empty(end)
    smoothed=np.append(add_start,smoothed)
    smoothed=np.append(smoothed,add_end)
    return smoothed


def smooth_butterworth(data, cutoff_freq, sampling_rate, order=4):
    nyquist = 0.5 * sampling_rate
    normal_cutoff = cutoff_freq / nyquist
    b, a = butter(order, normal_cutoff, btype="low", analog=False)
    y_smoothed = filtfilt(b, a, data)
    return y_smoothed