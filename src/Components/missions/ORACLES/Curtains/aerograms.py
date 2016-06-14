import os
import numpy as np
import datetime as dt
# matplotlib - MPLCONFIGDIR issue
try:
    import matplotlib as mpl
except:
    import tempfile
    import atexit
    import shutil

    mpldir = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, mpldir) #rm dir on exit

    os.environ['MPLCONFIGDIR'] = mpldir

    import matplotlib as mpl
# matplotlib display backend
mpl.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
# netCDF4 - PYTHON_EGG_CACHE issue
try:
    from netCDF4 import Dataset
except:
    import tempfile
    import atexit
    import shutil

    eggdir = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, eggdir) # rm dir on exit

    os.environ['PYTHON_EGG_CACHE'] = eggdir

    from netCDF4 import Dataset
import multiprocessing as mp

import find

#==============================================================================

cmap = {
    'meteo':plt.cm.Greens,
    'oc':plt.cm.bone_r,
    'bc':plt.cm.bone_r,
    'ss':plt.cm.Purples,
    'du':plt.cm.Reds,
    'su':plt.cm.pink_r,
    'co':plt.cm.YlGn,
    'all':plt.cm.copper_r,
}

width = 0.119
linewidth = 0.
global forecast

korus_cities = {
    '38.9x-77.0':'Washington_DC',
    '26.9x128.2':'citya',
    '32.1x125.5':'cityb',
    '32.7x128.7':'cityc',
    '33.2x126.2':'cityd',
    '33.2x130.3':'citye',
    '33.5x130.4':'cityf',
    '33.9x124.6':'cityg',
    '35.0x126.4':'cityh',
    '35.2x126.8':'cityi',
    '35.2x129.0':'cityj',
    '35.6x129.2':'cityk',
    '35.8x128.6':'cityl',
    '35.9x126.7':'citym',
    '36.0x127.0':'cityn',
    '36.3x127.4':'cityo',
    '36.5x126.3':'cityp',
    '37.0x127.0':'cityq',
    '37.3x127.3':'cityr',
    '37.4x126.9':'citys',
    '37.4x137.4':'cityt',
    '37.5x126.6':'cityu',
    '37.5x126.9':'cityv',
    '37.5x127.1':'cityw',
    '37.7x128.7':'cityx',
    '37.8x128.8':'cityy',
    '38.0x124.6':'cityz',
    '43.4x143.8':'cityzz'
}

#==============================================================================

def rc():
    '''custom plot settings, use mpl.rcdefaults() to default'''
    mpl.rc('lines', linewidth=0.5, antialiased=True)
    mpl.rc('patch', linewidth=0.5, facecolor='348ABD', edgecolor='eeeeee', \
          antialiased=True)
    mpl.rc('axes', facecolor='w', edgecolor='black', linewidth=0.5)
    mpl.rc('font', family='sans-serif', size=10.0)
    mpl.rc('xtick', color='black')
    mpl.rc('xtick.major', size=4, pad=6)
    mpl.rc('xtick.minor', size=2, pad=6)
    mpl.rc('ytick', color='black')
    mpl.rc('ytick.major', size=4, pad=6)
    mpl.rc('ytick.minor', size=2, pad=6)
    mpl.rc('legend', fancybox=True, fontsize=10.0)
    mpl.rc('figure', figsize='9.6, 9', dpi=100, facecolor='white') #88
    mpl.rc('figure.subplot', hspace=0.5, left=0.07, right=0.95, bottom=0.1, \
          top=0.95)

def custom_date_formatter(x, p):
    '''xaxis formatter for HHz + date + Year'''
    global forecast
    ti = mdates.num2date(x)
    if ti.hour == 0:
        if x < mdates.date2num(forecast + dt.timedelta(days=1)):
            return ti.strftime('%Hz\n%a %-d %b\n%Y')
        return ti.strftime('%Hz\n%a %-d %b')
    else:
        ti.strftime('%Hz')

def do_plots(ipath, stations, fcst, opath):
    '''receive dict of stations'''
    if not os.path.exists(ipath):
        return
    s = []
    for sta in stations:
        s.append((
            ipath, sta, float(stations[sta]['lat']),
            float(stations[sta]['lon']), fcst, opath
        ))

    #import time
    #start = time.time()
    # parallelize across stations
    pool = mp.Pool()
    pool.map(plot_wrapper, s)
    pool.close()
    pool.join()
    #print 'Elapsed Time: %.2fs' % ((time.time()-start))
    return

def plot_wrapper(args):
    return plot(*args)

def plot(ipath, station, lat, lon, fcst, opath):
    '''main plot driver'''
    global forecast
    forecast = fcst
    # plot settings
    rc()

    # read in all data for one station
    cldhgh = get_data(ipath, 'CLDHGH', station)
    cldmid = get_data(ipath, 'CLDMID', station)
    cldlow = get_data(ipath, 'CLDLOW', station)

    # could move to plotting routine since these depend upon product
    rh = get_data(ipath, 'RH', station)
    ocext = get_data(ipath, 'OCEXT', station)
    bcext = get_data(ipath, 'BCEXT', station)
    suext = get_data(ipath, 'SUEXT', station)
    co = get_data(ipath, 'CO', station)
    airdens = get_data(ipath, 'AIRDENS', station)
    ssext = get_data(ipath, 'SSEXT', station)
    duext = get_data(ipath, 'DUEXT', station)

    u = get_data(ipath, 'U', station)
    v = get_data(ipath, 'V', station)

    slp = get_data(ipath, 'SLP', station)
    t2m = get_data(ipath, 'T2M', station)

    ssexttau = get_data(ipath, 'SSEXTTAU', station)
    duexttau = get_data(ipath, 'DUEXTTAU', station)
    bcexttau = get_data(ipath, 'BCEXTTAU', station)
    ocexttau = get_data(ipath, 'OCEXTTAU', station)
    suexttau = get_data(ipath, 'SUEXTTAU', station)

    prectot = get_data(ipath, 'PRECTOT', station)
    precsno = get_data(ipath, 'PRECSNO', station)
    preccon = get_data(ipath, 'PRECCON', station)
    u2m = get_data(ipath, 'U2M', station)
    v2m = get_data(ipath, 'V2M', station)

    for product in ['meteo', 'oc', 'bc', 'ss', 'du', 'su', 'all', 'co']:
        # start plot template
        fig = plt.figure()
        ax1 = plt.subplot2grid((30, 1), (1, 0))
        ax2 = plt.subplot2grid((30, 1), (2, 0), sharex=ax1)
        ax3 = plt.subplot2grid((30, 1), (3, 0), sharex=ax1)
        ax4 = plt.subplot2grid((30, 1), (4, 0), rowspan=15, sharex=ax1)
        ax5 = plt.subplot2grid((30, 1), (21, 0), rowspan=4)
        ax7 = plt.subplot2grid((30, 1), (26, 0), rowspan=4)

        # CLDHGH
        cldhgh_t = [fcst+dt.timedelta(seconds=int(x)) for x in cldhgh['time']]
        ax1.set_axis_bgcolor('b')
        ax1.bar(cldhgh_t, np.absolute(cldhgh['data']), width=width, color='w', align='center', linewidth=linewidth)
        ax1.bar(cldhgh_t, np.negative(np.absolute(cldhgh['data'])), width=width, color='w', align='center', linewidth=linewidth)
        ax1.axhline(0, color='b')
        ax1.yaxis.set_ticks([])
        ax1.set_xlim([fcst, cldhgh_t[-1]])
        ax1.set_ylim(-1.5,1.5)
        ax1.xaxis.set_ticks_position('top')
        ax1.tick_params(axis=u'both', which=u'both', length=0)
        ax1.set_ylabel('High', rotation=0)
        ax1.yaxis.set_label_coords(-0.03, 0.1)
        del cldhgh_t

        # CLDMID
        cldmid_t = [fcst+dt.timedelta(seconds=int(x)) for x in cldmid['time']]
        ax2.set_axis_bgcolor('b')
        ax2.bar(cldmid_t, np.absolute(cldmid['data']), width=width, color='w', align='center', linewidth=linewidth)
        ax2.bar(cldmid_t, np.negative(np.absolute(cldmid['data'])), width=width, color='w', align='center', linewidth=linewidth)
        ax2.axhline(0, color='b')
        ax2.set_ylim(-2,2)
        ax2.yaxis.set_ticks([])
        ax2.tick_params(axis=u'both', which=u'both', length=0)
        ax2.set_ylabel('Mid', rotation=0)
        ax2.yaxis.set_label_coords(-0.03, 0.1)
        del cldmid_t
    
        # CLDLOW
        cldlow_t = [fcst+dt.timedelta(seconds=int(x)) for x in cldlow['time']]
        ax3.set_axis_bgcolor('b')
        ax3.bar(cldlow_t, np.absolute(cldlow['data']), width=width, color='w', align='center', linewidth=linewidth)
        ax3.bar(cldlow_t, np.negative(np.absolute(cldlow['data'])), width=width, color='w', align='center', linewidth=linewidth)
        ax3.axhline(0, color='b')
        ax3.yaxis.set_ticks([])
        ax3.tick_params(axis=u'both', which=u'both', length=0)
        ax3.set_ylim(-2,2)
        ax3.set_ylabel('Low', rotation=0)
        ax3.yaxis.set_label_coords(-0.03, 0.1)
        del cldlow_t
    
        # skip ax4, ax5, ax6 till product loop (3d fields + slp/t2m)
        fig.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in fig.axes[1:3]+fig.axes[4:-1]], visible=False)
    
        # precip
        precip_t = [fcst+dt.timedelta(seconds=int(x)) for x in prectot['time']]
        ax7.bar(precip_t, prectot['data']*10000, width=width, color='g', align='center', linewidth=linewidth)
        ax7.bar(precip_t, precsno['data']*10000, width=width*.6, color='b', align='center', linewidth=linewidth)
        ax7.bar(precip_t, preccon['data']*10000, width=width*.4, color='r', align='center', linewidth=linewidth)
        ylow, yhigh = ax7.get_ylim()
        ax7.locator_params(axis='y', nbins=3)
        ax7.yaxis.grid(True)
    
        # 2m winds
        ax8 = ax7.twinx()
        _2m_t = [fcst+dt.timedelta(seconds=int(x)) for x in u2m['time']]
        _2m = np.sqrt(u2m['data']*u2m['data'] + v2m['data']*v2m['data'])*1.943846
        ax8.plot(_2m_t, _2m, color='k')
        ax8.set_xlim([fcst, _2m_t[-1]])
        ax8.locator_params(axis='y', nbins=3)
        ax7.xaxis.set_ticks([x for x in precip_t if x.hour == 0 or x.hour == 12])
        ax7.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(custom_date_formatter))
        del _2m_t
    
        # ax4, ax5, ax6 here
        rh_t = [fcst+dt.timedelta(seconds=int(t)) for t in rh['time']]
        ax4.xaxis.set_ticks([t for t in rh_t if t.hour == 0 or t.hour == 12])
        ax4.xaxis.set_major_formatter(mdates.DateFormatter('%Hz'))

        x = np.arange(len(rh_t))
        x = [fcst+dt.timedelta(hours=int(ti) * 3) for ti in x] # 3-hourly data!
        x = mdates.date2num(x)
        xi = x

        # 3d field
        if product in 'meteo':
            y = [i for i in range(len(rh['lev'])) if 600. <= rh['lev'][i] <= 950.]
            ylabels = [int(rh['lev'][i]) for i in y]
            yi = y
            x, y = np.meshgrid(xi, y)
            rh_v = rh['data'].transpose()[yi]*100
            ax4.contourf(x, y, rh_v, cmap=cmap[product])
            CS = ax4.contour(x, y, rh_v, colors='k', levels=[20,40,60,70,80,85,90,95])
            fmt = '%d'
        else:
            if product in 'all':
                y = [i for i in range(len(ocext['lev'])) if 300. <= ocext['lev'][i] <= 950.]
                ylabels = [int(ocext['lev'][i]) for i in y]
                yi = y
                x, y = np.meshgrid(xi, y)
                _3d = ocext['data'].transpose()[yi] + bcext['data'].transpose()[yi] + suext['data'].transpose()[yi]
                _3d *= 1000 # m-1 to km-1
                ax4.contourf(x, y, _3d, cmap=cmap[product])
                CS = ax4.contour(x, y, _3d, colors='k')
                fmt = '%.2f'
            elif product in 'co':
                y = [i for i in range(len(co['lev'])) if 300. <= co['lev'][i] <= 950.]
                ylabels = [int(co['lev'][i]) for i in y]
                yi = y
                x, y = np.meshgrid(xi, y)
                _3d = co['data'].transpose()[yi] * 1e7 * airdens['data'].transpose()[yi]
                ax4.contourf(x, y, _3d, cmap=cmap[product])
                CS = ax4.contour(x, y, _3d, colors='k')
                fmt = '%.2f'
            else:
                _3d = {
                    'oc':ocext, 'bc':bcext, 'ss':ssext, 'du':duext, 'su':suext
                }
                y = [i for i in range(len(_3d[product]['lev'])) if 300. <= _3d[product]['lev'][i] <= 950.]
                ylabels = [int(_3d[product]['lev'][i]) for i in y]
                yi = y
                x, y = np.meshgrid(xi, y)
                ax4.contourf(x, y, _3d[product]['data'].transpose()[yi]*1000, cmap=cmap[product])
                CS = ax4.contour(x, y, _3d[product]['data'].transpose()[yi]*1000, colors='k')
                if product in ['oc', 'du', 'su']:
                    fmt = '%.3f'
                else:
                    fmt = '%.4f'
        plt.setp(CS.collections, linewidth=0.2)
        plt.setp(CS.collections, linestyle='--')
        tl = ax4.clabel(CS, CS.levels, inline=1, fontsize=8, fmt=fmt)
        for te in tl:
            te.set_bbox(dict(color='w', alpha=0.8, pad=0.1))
        ax4.set_yticklabels(ylabels[::2][:-1])
        ax4.yaxis.set_tick_params(direction='out')
        ax4.yaxis.tick_left()
        ax4.xaxis.set_tick_params(direction='out')
        ax4.xaxis.tick_bottom()
        ax4.set_ylim(yi[1], yi[-1])

        # u/v winds
        ax4.barbs(
            x, y, u['data'].transpose()[yi]*2.23694, v['data'].transpose()[yi]*2.23694,
            barb_increments=dict(half=5, full=10, flag=50)
        )

        # slp & t2m
        if product in 'meteo':
            slp_t = [fcst+dt.timedelta(seconds=int(x)) for x in slp['time']]
            ax5.plot(slp_t, slp['data']/100, color='k')
            ax5.yaxis.grid(True)
            ax6 = ax5.twinx()
            ax6.plot(slp_t, (t2m['data']-273.15)*1.8+32., color='r')
            ax5.tick_params(axis=u'both', which=u'both', length=0)
            ax5.locator_params(axis='y', nbins=3)
            ax6.locator_params(axis='y', nbins=3)
            for t1 in ax6.get_yticklabels():
                t1.set_color('r')
        else:
            tau_t = [fcst+dt.timedelta(seconds=int(x)) for x in suexttau['time']]
            _su = suexttau['data']
            _c = _su + bcexttau['data'] + ocexttau['data']
            _du = _c + duexttau['data']
            _ss = _du + ssexttau['data']
            ax5.bar(tau_t, _ss, color='b', width=width*.2, align='center', linewidth=linewidth)
            ax5.bar(tau_t, _du, color='r', width=width*.2, align='center', linewidth=linewidth)
            ax5.bar(tau_t, _c, color='k', width=width*.2, align='center', linewidth=linewidth)
            ax5.bar(tau_t, _su, color='y', width=width*.2, align='center', linewidth=linewidth)
            ax5.yaxis.grid(True)
            ax5.tick_params(axis=u'both', which=u'both', length=0)
            ax5.locator_params(axis='y', nbins=3)

        # Text
        te = fig.text(-0.03, 0.88, 'Clouds (%)', va='center', rotation='vertical', fontsize=11, color='w')
        te.set_bbox(dict(color='b', alpha=1))
        te = fig.text(0., 0.6, 'Wind Barbs (mi/hr)', va='center', rotation='vertical', fontsize=11)
        te_ = fig.text(-0.06, 0.15, 'Rain', va='center', rotation='vertical', fontsize=11)
        te_.set_bbox(dict(color='lawngreen', alpha=1))
        te = fig.text(-0.03, 0.15, 'Snow (Liq)', va='center', rotation='vertical', fontsize=11, color='w')
        te.set_bbox(dict(color='b', alpha=1))
        te = fig.text(0., 0.15, 'Convective', va='center', rotation='vertical', fontsize=11)
        te.set_bbox(dict(color='r', alpha=1))
        te = fig.text(0.03, 0.15, '(mm)', va='center', rotation='vertical', fontsize=11)
        te_1 = fig.text(1., 0.15, '2m w spd', va='center', rotation=-90, fontsize=11)
        te = fig.text(0.98, 0.15, '(knots)', va='center', rotation=-90, fontsize=11)

        if product in 'meteo':
            te = fig.text(-0.03, 0.6, 'Relative Humidity (%)', va='center', rotation='vertical', fontsize=11)
            te.set_bbox(dict(color='lawngreen', alpha=1))
            te = fig.text(0., 0.3, 'T2M (F)', va='center', rotation='vertical', fontsize=11, color='r')
            te = fig.text(-0.03, 0.3, 'SLP (hPa)', va='center', rotation='vertical', fontsize=11)
        else:
            if product in 'su':
                te = fig.text(-0.03, 0.6, 'SO$_4$ Extinction (1/km)', va='center', rotation='vertical', fontsize=11)
            elif product in 'all':
                te = fig.text(-0.06, 0.6, 'Org Carbon+Blk Carbon+SO$_4$\nAerosol Extinction (1/km)', va='center', rotation='vertical', fontsize=11)
            elif product in 'co':
                te = fig.text(-0.03, 0.6, 'CO concentration (x10$^7$)', va='center', rotation='vertical', fontsize=11)
            else:
                te = fig.text(-0.03, 0.6, product.upper()+' Extinction (1/km)', va='center', rotation='vertical', fontsize=11)
            te.set_bbox(dict(color=mpl.cm.get_cmap(cmap[product])(0.35), alpha=1))
            te = fig.text(-0.04, 0.3, 'Sea Salt', va='center', rotation='vertical', fontsize=11, color='b')
            te = fig.text(-0.02, 0.3, 'Dust', va='center', rotation='vertical', fontsize=11, color='r')
            te = fig.text(0., 0.3, 'Carbon', va='center', rotation='vertical', fontsize=11, color='k')
            te = fig.text(0.02, 0.3, 'SO4', va='center', rotation='vertical', fontsize=11, color='y')
        lbl = fig.text(0.5,0.01, 'Lat = %s, Lon = %s, Location = %s, Fcst_Init = %s' % (lat, lon, station, forecast), ha='center', fontsize=12)

        #img = '.'.join(['forecast',product,fcst.strftime('%Y%m%d'), str(lat), str(lon)])+'.png'
        img = '.'.join(['forecast', 'meteogram', forecast.strftime('%Y%m%d'), str(lat), str(lon), korus_cities[str(lat)+'x'+str(lon)]])+'.png'
        plt.savefig(os.path.join(opath, img), bbox_inches='tight', dpi=100, bbox_extra_artists=(lbl,te_,te_1))

        plt.close()
    del cldhgh, cldmid, cldlow, u, v, slp, t2m
    del rh, ocext, bcext, suext, co, airdens, ssext, duext
    del ssexttau, duexttau, bcexttau, ocexttau, suexttau
    del prectot, precsno, preccon, u2m, v2m

def get_data(ipath, field, station):
    try:
        fi = [x for x in find.find(path=ipath, ext='.nc') if '_'+field+'.' in x]
        if fi:
            fi = fi[0]
            d = Dataset(fi, 'r')

            # get station index
            st = [(i,x) for i,x in enumerate(d.ncattrs()) if 'Station' in x]
            s = [x for i,x in st if getattr(d, x) in station][0]
            s = int(s.split('_')[-1]) - 1

            if 'time' not in d.variables:
                return None
            # have to slice instead of obtaining pointer to netCDF4.Variable
            time = d.variables['time'][:]

            if 'lev' in d.variables:
                lev = d.variables['lev'][:]
            else:
                lev = None

            if field not in d.variables:
                return None
            data = d.variables[field][s]
            d.close()
            del d, s, st, fi
            return {'name':field, 'data':data, 'time':time, 'lev':lev}
    except:
        return None
