import numpy as np
from scipy.stats import gmean
from sklearn.preprocessing import StandardScaler

def thin_layers_removal_old(h_min, md_list):
    """
    Removes thin layers from the input list based on the specified minimum height difference.
    
    Parameters:
        h_min: The minimum height difference required for layer removal.
        md_list: The list of heights to be processed.
        
    Returns:
        The modified list after removing the thin layers.
    """
    ZN_list = []
    if h_min > 0:
        i = 0
        while i < len(md_list) - 2:
            print(md_list[i])
            if (md_list[i+1] - md_list[i]) < h_min:
                if (md_list[i+1] - md_list[i]) > (md_list[i+2] - md_list[i+1]) and ((md_list[i+1] - md_list[i]) + (md_list[i+2] - md_list[i+1])) >= h_min:
                    del md_list[i+1]
                    if i+1 < len(ZN_list):  # Check if index is within ZN_list length
                        del ZN_list[i+1]
                elif (md_list[i] - md_list[i-1]) >= (md_list[i+2] - md_list[i+1]):
                    del md_list[i]
                    if i < len(ZN_list):  # Check if index is within ZN_list length
                        del ZN_list[i]
                else:
                    del md_list[i+1]
                    if i+1 < len(ZN_list):  # Check if index is within ZN_list length
                        del ZN_list[i+1]
            else:
                i += 1
    return md_list


import numpy as np

def thin_layers_removal(h_min, md_list):
    """
    Removes thin layers from the input list based on the specified minimum height difference.
    
    Parameters:
        h_min: The minimum height difference required for layer removal.
        md_list: The list of heights to be processed.
        
    Returns:
        The modified list after removing the thin layers.
    """
    if h_min > 0:
        md_list = list(md_list)  # Ensure md_list is a list
        i = 0
        while i < len(md_list) - 1:
            if (md_list[i+1] - md_list[i]) < h_min:
                if i + 2 < len(md_list) and (md_list[i+1] - md_list[i]) > (md_list[i+2] - md_list[i+1]) and ((md_list[i+1] - md_list[i]) + (md_list[i+2] - md_list[i+1])) >= h_min:
                    del md_list[i+1]
                elif i > 0 and (md_list[i] - md_list[i-1]) >= (md_list[i+2] - md_list[i+1]):
                    del md_list[i]
                else:
                    del md_list[i+1]
            else:
                i += 1
        md_list = np.array(md_list)  # Convert list back to numpy array
    return md_list



def layering(variable, H_min, percentile, sensitivity):
    """
    Computes layers based on the input variable, H_min, percentile, and sensitivity. 
    Parameters:
        variable: The input variable containing reference data.
        H_min: The minimum height difference required for layering.
        percentile: The percentile value for computing dmax and dmin.
        sensitivity: The sensitivity factor for computing sens_max and sens_min.
    Returns:
        List of depths representing the layered structure.
    """
    asc_list = []
    des_list = []
    A_list = []
    B_list = []
    ZN_list = []
    
    asc_flag = 0
    des_flag = 0
    
    lc = variable.reference.size
    md = variable.reference

    sr = (variable.reference.max()-variable.reference.min())/variable.reference.size
    der = np.gradient(variable.first_data_column, round(sr, 2))

    dmax = np.nanpercentile(der, 100-percentile)
    dmin = np.nanpercentile(der, percentile)

    sens_max = dmax*(1-sensitivity)
    sens_min = abs(dmin*(1-sensitivity))
    
    for i in range(lc):
        if der[i] > 0 and der[i] >= sens_max:
            asc_list.append(der[i])
            if asc_flag == 0:
                i_start = i
                asc_flag = 1
        else:
            if asc_flag == 1:
                if len(A_list) == 0 or (md[int(i_start + (i - i_start) // 2)] - A_list[-1]) >= H_min:
                    A_list.append(md[int(i_start + (i - i_start) // 2)])
                    ZN_list.append(1)
                asc_list = []
                asc_flag = 0

    for i in range(lc):
            if der[i] < 0 and abs(der[i]) >= sens_min:
                des_list.append(abs(der[i]))
                if des_flag == 0:
                    i_start = i
                    des_flag = 1
            else:
                if des_flag == 1:
                    if len(B_list) == 0 or (md[int(i_start + (i - i_start) // 2)] - B_list[-1]) >= H_min:
                        B_list.append(md[int(i_start + (i - i_start) // 2)])
                        ZN_list.append(1)
                    des_list = []
                    des_flag = 0

    md_list = A_list + B_list
    md_list = sorted(md_list)
    
    return md_list



def calculate_averages(MD_list, variable, averaging):
    """
    Calculates averages based on the input parameters MD_list, variable, and averaging.
    
    Parameters:
        MD_list: List of measured depths to calculate averages.
        variable: The input variable containing reference data.
        averaging: The type of averaging to perform (e.g., Arithmetic, Geometric, Median, Realistic).
    
    Returns:
        numpy array: An array containing the calculated averages.
    """
    averages = []
    gr_val = variable.first_data_column
    md = variable.reference
    lc = variable.reference.size
    if averaging == 'None':
        return np.array([0])
    
    z = 0
    while z < len(MD_list)-1:
        interval_values = []
        
        for i in range(lc):
            if gr_val[i] is not None:
                if MD_list[z] <= md[i] < MD_list[z+1]:
                    interval_values.append(gr_val[i])
                
                if md[i] >= MD_list[z+1] or i == lc - 1:
                    if averaging == 'Arithmetic':
                        averages.append(np.mean(interval_values))
                    elif averaging == 'Geometric':
                        averages.append(gmean(interval_values))
                    elif averaging == 'Frequency':
                        averages.append(np.argmax(np.bincount(interval_values)))
                    elif averaging == 'Median':
                        averages.append(np.median(interval_values))
                    elif averaging == 'Realistic':
                        INFL = 0
                        for x in range(len(interval_values)):
                            if x != 0 and x != len(interval_values)-1:
                                if interval_values[x] < interval_values[x-1] and interval_values[x+1] >= interval_values[x]:
                                    if INFL != 2:
                                        INFL = 1  # low extremum
                                    else:
                                        INFL = 0
                                        break
                                if interval_values[x] >= interval_values[x-1] and interval_values[x+1] < interval_values[x]:
                                    if INFL != 1:
                                        INFL = 2  # high extremum
                                    else:
                                        INFL = 0
                                        break
                        if INFL == 1:
                            averages.append(min(interval_values))
                        if INFL == 2:
                            averages.append(max(interval_values))
                        if INFL == 0:
                            averages.append(np.mean(interval_values))
                    z += 1
                    interval_values = []
                    break

    
    return np.append(np.array(averages), 0)

def combine_in_one_log(inter_logs, option=1):
    """
    Combine multiple logs into a single log by stacking them horizontally. 
    If option is 1, logs are stacked directly, otherwise, the absolute values are stacked. 
    Calculate the common log by taking the mean along the columns. 
    Return the common log.
    """
    scaler = StandardScaler()
    if option == 1:
        all_logs = np.column_stack([scaler.fit_transform(log.values) for log in inter_logs])
    else:
        all_logs = np.column_stack([abs(scaler.fit_transform(log.values)) for log in inter_logs])
    common_log = np.nanmean(all_logs, axis=1)
    return common_log
