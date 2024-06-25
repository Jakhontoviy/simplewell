import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
import plotly.express as px

def readXLSX(filename):
    #file_path = os.path.join(os.path.dirname(__file__), 'source', filename)
    df = pd.read_excel(filename, header=0)
    return df

def chart1d_model(log, path):
    chart = readXLSX(path)
    model = np.poly1d(np.polyfit(chart['X'], chart['Y'], 4))
    result = model(log)
    return result

def chart2d_model_v1(x_log_ini, y_log_ini, path, show_chart=False):
    chart = readXLSX(path)
    x_chart = np.interp(chart['X'], (chart['X'].min(), chart['X'].max()), (-1, +1))
    y_chart = np.interp(chart['Y'], (chart['Y'].min(), chart['Y'].max()), (-1, +1))
    x_log = np.interp(x_log_ini, (chart['X'].min(), chart['X'].max()), (-1, +1))
    y_log = np.interp(y_log_ini, (chart['Y'].min(), chart['Y'].max()), (-1, +1))
    z_chart = chart['Z']
    rbf = Rbf(x_chart, y_chart, z_chart, function="linear")
    z_log_pred = rbf(x_log, y_log)
    if show_chart == True:
        log_points = {
            'X': x_log_ini,
            'Y': y_log_ini,
            'Z': z_log_pred
        }
        plot_2d_chart(chart, log_points)
    return z_log_pred

def chart2d_model(x_log_ini, y_log_ini, path, show_chart=False):
    chart = readXLSX(path)
    x_min, x_max = chart['X'].min(), chart['X'].max()
    y_min, y_max = chart['Y'].min(), chart['Y'].max()
    x_chart = rescale_linear(chart['X'], -1, +1, x_min, x_max)
    y_chart = rescale_linear(chart['Y'], -1, +1, y_min, y_max)
    x_log = rescale_linear(x_log_ini, -1, +1, x_min, x_max)
    y_log = rescale_linear(y_log_ini, -1, +1, y_min, y_max)
    z_chart = chart['Z']
    rbf = Rbf(x_chart, y_chart, z_chart, function="linear")
    z_log_pred = rbf(x_log, y_log)
    if show_chart == True:
        log_points = {
            'X': x_log_ini,
            'Y': y_log_ini,
            'Z': z_log_pred
        }
        plot_2d_chart(chart, log_points)
    return z_log_pred

def plot_2d_chart(chart_data, log_data):
    fig = px.line(chart_data, x=chart_data['X'], y=chart_data['Y'], color = chart_data['Z'], markers = True)
    fig.add_scatter(x=log_data['X'], y=log_data['Y'], mode='markers', marker=dict(color=log_data['Z']), text=log_data['Z'])
    fig.show()


def rescale_linear(array, new_min, new_max, minimum, maximum):
    """Rescale an arrary linearly."""
    #minimum, maximum = np.min(array), np.max(array)
    m = (new_max - new_min) / (maximum - minimum)
    b = new_min - m * minimum
    return m * array + b