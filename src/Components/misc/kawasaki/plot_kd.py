"""
Plot KD time series.
"""

from matplotlib.pyplot import *
from numpy import *
from npz   import NPZ

def _P(x):
    return 100.*x/x.sum()

def plot_decades():

    figure()

    n = NPZ('NPZ/kawasaki_japan.monthly.1969-2010.npz')
    Year  = array([t.year  for t in n.tyme])
    Month = array([t.month for t in n.tyme])

    I79 = (Year==1979)
    I82 = (Year==1982)
    I86 = (Year==1986)
    
    c70 = NPZ('NPZ/kawasaki_japan.climatology.1969-1979.npz')
    c80 = NPZ('NPZ/kawasaki_japan.climatology.1980-1989.npz')
    c90 = NPZ('NPZ/kawasaki_japan.climatology.1990-1999.npz')
    c00 = NPZ('NPZ/kawasaki_japan.climatology.2000-2010.npz')

    month = list(range(12))
    figure(dpi=120)

    # Full counts
    # -----------
    plot(month,c00.all,'-o',linewidth=1.5,label='2000-2010')
    plot(month,c90.all,'-o',linewidth=1.5,label='1990-1999')
    plot(month,c80.all,'-o',linewidth=1.5,label='1980-1989')
    plot(month,c70.all,'-o',linewidth=1.5,label='1969-1979')

    xticks(month,['Jan','Feb','Mar','Apr','May','Jun','Jul',
                  'Aug','Sep','Oct','Nov','Dec'])
    ylabel('No. of Cases/Month')
    #xlabel('Climatological Month')
    title('Climatological Number of KD Cases')
    grid()
    axis([-0.5,11.5,0,1300])
    legend(loc='upper center',prop={'size':10})

    ax2=twinx()
    ax2.plot(month,n.all[I79],'k:o',linewidth=1.5,label='1979')
    ax2.plot(month,n.all[I82],'y:o',linewidth=1.5,label='1982')
    ax2.plot(month,n.all[I86],'m:o',linewidth=1.5,label='1986')
    legend(loc='upper right',prop={'size':10})
    axis([-0.5,11.5,0,3500])

    savefig('Images/kd.climatology.counts.png',
            bbox_inches='tight',dpi=120)

    # % counts
    # -----------
    figure(dpi=120)
    plot(month,_P(c00.all),'-o',linewidth=1.5,label='2000-2010')
    plot(month,_P(c90.all),'-o',linewidth=1.5,label='1990-1999')
    plot(month,_P(c80.all),'-o',linewidth=1.5,label='1980-1989')
    plot(month,_P(c70.all),'-o',linewidth=1.5,label='1969-1979')

    xticks(month,['Jan','Feb','Mar','Apr','May','Jun','Jul',
                  'Aug','Sep','Oct','Nov','Dec'])
    ylabel('%')
    #xlabel('Climatological Month')
    title('Climatological % Number of KD Cases per Month')
    grid()
    axis([-0.5,11.5,0,13])
    legend(loc='lower left',prop={'size':10})

    ax2=twinx()
    ax2.plot(month,_P(n.all[I79]),'k:o',linewidth=1.5,label='1979')
    ax2.plot(month,_P(n.all[I82]),'y:o',linewidth=1.5,label='1982')
    ax2.plot(month,_P(n.all[I86]),'m:o',linewidth=1.5,label='1986')
    axis([-0.5,11.5,0,25])
    legend(loc='upper center',prop={'size':10})
    
    savefig('Images/kd.climatology.percent.png',
            bbox_inches='tight',dpi=120)

def plot_ts():

    d = NPZ('NPZ/kawasaki_japan.daily.1969-2010.npz')
    m = NPZ('NPZ/kawasaki_japan.monthly.1969-2010.npz')
    y = NPZ('NPZ/kawasaki_japan.yearly.01-12.npz')

    figure(dpi=120)
    subplot(211)
    plot(d.tyme,d.all,label='daily')
    grid()
    title('Kawasaki Disease Cases')
    ylabel('Total No. per Day')
    
    subplot(212)
    plot(m.tyme,m.all,'k',label='per month')
    plot(y.tyme,y.all,'r',linewidth=2,label='per year')
    ylabel('Average Number')
    grid()
    legend(loc='upper right',prop={'size':10})
    
    savefig('Images/kd.timeseries.png',
            bbox_inches='tight',dpi=120)

def plot_seasonal():
    from mk_npz import yearly
    
    figure(dpi=120,figsize=(11,11*9./16))
    
    tall, all = yearly(m1=1,m2=12,npzWrite=False)
    tfma, fma = yearly(m1=2,m2=4,npzWrite=False)
    tmjj, mjj = yearly(m1=5,m2=7,npzWrite=False)
    taso, aso = yearly(m1=8,m2=10,npzWrite=False)
    tndj, ndj = yearly(m1=11,m2=13,npzWrite=False)

    year = list(range(1969,2011))
    
    plot(year,ndj,'-',linewidth=1.5,label='NOV-JAN')
    plot(year,fma,'-',linewidth=1.5,label='FEB-MAY')
    plot(year,mjj,'-',linewidth=1.5,label='MAY-JUL')
    plot(year,aso,'-',linewidth=1.5,label='AUG-OCT')
    plot(year,all,'k-',linewidth=3,label='YEAR')
    grid()
    legend(loc='upper right',prop={'size':10})
    title('Kawasaki Disease Cases: Seasons by Year')
    
    savefig('Images/kd.seasonal.png',
            bbox_inches='tight',dpi=120)

def plot_population():

    # Japanese Population, 1969-2009
    # Source: http://www.stat.go.jp/english/data/chouki/02.htm
    # --------------------------------------------------------
    Pop04 = array([1573824,1898242,1883710,1839328,1843664,1442620,1913757,2028786,2063827,2020504,1973733,1586504,1636235,1706055,1749163,1837459,1429658,1488363,1504910,1508715,1527617,1213685,1260478,1301517,1343438,1373779,1191578,1201008,1185263,1208065,1209340,1171652,1166160,1192157,1189303,1184826,1056800,1091316,1115649,1149450,1164872,NaN],dtype=float)

    y = NPZ('NPZ/kawasaki_japan.yearly.01-12.npz')
    year = list(range(1969,2011))

    figure(dpi=120,figsize=(11,11*9./16))
    plot(year,y.all/1000.,'b-o',linewidth=1.5,label='KD Cases (K)')
    plot(year,Pop04/1000000.,'r-o',linewidth=1.5,label='0-4 Years Old Population (M)')
    grid()
    legend(loc='upper right',prop={'size':10})
    title('Kawasaki Disease Cases by Year')
    
    savefig('Images/kd.population.png',
            bbox_inches='tight',dpi=120)
    

#....................................................................

if __name__ == "__main__":

    plot_population()
    
