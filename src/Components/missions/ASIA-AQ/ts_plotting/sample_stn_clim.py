
from pyobs.sampler import STATION
from glob import glob

Variables = ['BCEXTTAU','BCSMASS','DUEXTTAU','DUSMASS','DUSMASS25',
             'OCEXTTAU','OCSMASS','SO4SMASS','SSEXTTAU','SSSMASS',
             'SSSMASS25','SUEXTTAU','TOTEXTTAU']

if __name__ == "__main__":

    lons = [126.99]
    lats = [37.55]
    stations = 'seoul'

#    lons = [120.98]
#    lats = [14.599]
#    stations = 'manila'

    years = '2014',
#    years = '2015','2016','2017','2018','2019','2020','2021','2022','2023'
#    years = '2020','2021','2022','2023'

    for year in years:
        dataset = sorted(glob('MERRA2_all/Y{}/M01/MERRA2.tavg1_2d_aer_Nx.{}*'.format(year,year)))
        dataset += sorted(glob('MERRA2_all/Y{}/M02/MERRA2.tavg1_2d_aer_Nx.{}*'.format(year,year)))
        stn = STATION([stations],lons,lats,dataset,verbose=True)

        ds = stn.sample(Variables=Variables)
        ds.to_netcdf(path='{}_{}.nc'.format(stations,year))
