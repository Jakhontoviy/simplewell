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