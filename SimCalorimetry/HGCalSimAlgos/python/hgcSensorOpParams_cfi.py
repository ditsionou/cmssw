import FWCore.ParameterSet.Config as cms

def hgcSiSensorIleak(version):

    """ 
    this method returns different parameterizations of the leakage current for different versions
    {600V,800V}_annealing - 2023 parameterizations for operations at -30C, 600V or 800V and different annealing times
    if version is unknown a ValueError exception is raised
    """
    
    if version=='600V_90m':
        return [1.0,-42.715]
    elif version=='800V_90m':
        return [1.0,-42.509]
    
    raise ValueError('version={} is unknown to retrieve Ileak parameterization for HGC Si sensors'.format(version))


def hgcSiSensorCCE(sensor,version):

    """ 
    this method returns different parameterizations of the charge collection efficiency (CCE)
    for different sensor types (sensor) and measurement versions (version)
    sensor = 120,200,300
    version = {600V,800V}_{nom,up,dn}_{annealing}   - 2023 base measurements at different voltages
    if the pair (sensor,version) is unknown a ValueError exception is raised
    """
    
    if version=='600V_nom_90m':
        if sensor==120  : return [-29.96445011, 1162.10274299]
        elif sensor==200: return [-32.60721635, 1216.43542273]
        elif sensor==300: return [-22.24623809, 824.65456822]
    if version=='600V_up_90m':
        if sensor==120  : return [-34.16741790170276, 1317.654494984719]
        elif sensor==200: return [-38.010147623701236, 1409.72800745149]
        elif sensor==300: return [-27.430094132194704, 1006.7703542871437]
    if version=='600V_dn_90m':
        if sensor==120  : return [-25.761482326365066, 1006.550990987932]
        elif sensor==200: return [-27.204285071882147, 1023.1428380033524]
        elif sensor==300: return [-17.06238204833479, 642.5387821618086]
    if version=='800V_nom_90m':
        if sensor==120  : return [-26.31895177, 1038.89425047]
        elif sensor==200: return [-27.35815924, 1039.49193707]
        elif sensor==300: return [-23.35206332, 876.10888287]
    if version=='800V_up_90m':
        if sensor==120  : return [-30.52199443454662, 1194.4487739797887]
        elif sensor==200: return [-33.68392908416032, 1265.1254176082737]
        elif sensor==300: return [-29.851729815283992, 1104.2151379728784]
    if version=='800V_dn_90m':
        if sensor==120  : return [-22.11590911057771, 883.3397269553066]
        elif sensor==200: return [-21.03238939261077, 813.8584565282233]
        elif sensor==300: return [-16.852396820599186, 648.0026277649256]
        

    raise ValueError('sensor={} version={} is unknown to retrieve CCE parameterization for HGC Si sensors'.format(sensor,version))
