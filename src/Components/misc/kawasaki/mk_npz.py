"""
Converts Kawasaki Disease data to npz format.
"""

from datetime import datetime, timedelta
from numpy import savez, array
from pyobs.npz import NPZ

def to_npz():
    
    Lines = Lines = open('KDjapan1969-2010dlyregct.txt').readlines()
    Header = Lines[0:4]
    Lines = Lines[4:]

    Tyme, Pref, All, Fem, Male = [], [], [], [], []
    
    for line in Lines:

        line = line.replace('\n','')

        y, m, d, pref, all, fem, male = line.split()

        Tyme.append(datetime(int(y),int(m),int(d)))
        Pref.append(int(pref))
        All.append(int(all))
        Fem.append(int(fem))
        Male.append(int(male))
        
    savez('NPZ/kawasaki_japan.obs.1969-2010.npz',
          Header=Header,
          tyme=array(Tyme),
          pref=array(Pref,dtype='int8'),
          all=array(All,dtype='int8'),
          fem=array(Fem,dtype='int8'),
          male=array(Male,dtype='int8'),
         )

def daily():

    day = timedelta(seconds=24*60*60)
    
    n = NPZ('NPZ/kawasaki_japan.obs.1969-2010.npz')

    t0, tf = n.tyme[0], n.tyme[-1]

    t = t0
    Tyme, All, Male, Fem = [], [], [], []
    while t < tf:

         I = n.tyme==t

         if any(I):
             all = n.all[I].sum()
             male = n.male[I].sum()
             fem = n.fem[I].sum()
         else:
             all, male, fem = 0, 0, 0

         Tyme.append(t)
         All.append(all)
         Male.append(male)
         Fem.append(fem)

         print('[] Got %s with %5d cases '%(str(t),all))
         
         t += day

    savez('kawasaki_japan.daily.1969-2010.npz',
          Header=n.Header,
          tyme=array(Tyme),
          all=array(All,dtype='int16'),
          fem=array(Fem,dtype='int16'),
          male=array(Male,dtype='int16'),
         )

def monthly():

    day = timedelta(seconds=24*60*60)
    
    n = NPZ('NPZ/kawasaki_japan.obs.1969-2010.npz')

    Year  = array([t.year  for t in n.tyme])
    Month = array([t.month for t in n.tyme])
    
    t0, tf = n.tyme[0], n.tyme[-1]

    t = t0
    Tyme, All, Male, Fem = [], [], [], []
    for year in range(t0.year,tf.year+1):
        for month in range(1,13):

            t = datetime(year,month,15)
            I = (Year==year)&(Month==month)

            if any(I):
                all = n.all[I].sum()
                male = n.male[I].sum()
                fem = n.fem[I].sum()
            else:
                all, male, fem = 0, 0, 0

            Tyme.append(t)
            All.append(all)
            Male.append(male)
            Fem.append(fem)

            print('[] Got %s with %5d cases '%(str(t),all))
         
    savez('NPZ/kawasaki_japan.monthly.1969-2010.npz',
          Header=n.Header,
          tyme=array(Tyme),
          all=array(All,dtype='int16'),
          fem=array(Fem,dtype='int16'),
          male=array(Male,dtype='int16'),
         )

def climatology(y1=1969,y2=2010):
    
    n = NPZ('NPZ/kawasaki_japan.monthly.1969-2010.npz')

    Year  = array([t.year  for t in n.tyme])
    Month = array([t.month for t in n.tyme])
    
    t0, tf = n.tyme[0], n.tyme[-1]

    t = t0
    Tyme, All, Male, Fem = [], [], [], []
    for month in range(1,13):

        t = datetime(1900,month,15)
        I = (Month==month)&(Year!=1979)&(Year!=1982)&(Year!=1986)\
                          &(Year>=y1)&(Year<=y2)

        if any(I):
            all = n.all[I].mean()
            male = n.male[I].mean()
            fem = n.fem[I].mean()
        else:
            all, male, fem = 0, 0, 0

        Tyme.append(t)
        All.append(all)
        Male.append(male)
        Fem.append(fem)

        print('[] Got %s with %5d cases '%(str(t),all+male+fem))
         
    savez('NPZ/kawasaki_japan.climatology.%d-%d.npz'%(y1,y2),
          Header=n.Header,
          tyme=array(Tyme),
          all=array(All,dtype='int16'),
          fem=array(Fem,dtype='int16'),
          male=array(Male,dtype='int16'),
         )

def yearly(m1=1,m2=12,npzWrite=True):
    
    n = NPZ('NPZ/kawasaki_japan.monthly.1969-2010.npz')

    Year  = array([t.year  for t in n.tyme])
    Month = array([t.month for t in n.tyme])

    if m2==13:
        Month[Month==1] = 13
    
    t0, tf = n.tyme[0], n.tyme[-1]
    
    t = t0
    Tyme, All, Male, Fem = [], [], [], []
    for year in range(t0.year,tf.year+1):

        halfy = (datetime(year,12,31) - datetime(year,1,1))/2 
        t = datetime(year,1,1) + halfy

        I = (Year==year)&(Month>=m1)&(Month<=m2)

        if any(I):
            all = n.all[I].mean()
            male = n.male[I].mean()
            fem = n.fem[I].mean()
        else:
            all, male, fem = 0, 0, 0

        Tyme.append(t)
        All.append(all)
        Male.append(male)
        Fem.append(fem)

        print('[] Got %s with %5d cases '%(str(t),all))

    tyme=array(Tyme),
    all=array(All,dtype='int16')
    fem=array(Fem,dtype='int16')
    male=array(Male,dtype='int16')

    if npzWrite:
        savez('NPZ/kawasaki_japan.yearly.%02d-%02d.npz'%(m1,m2),
              Header=n.Header,tyme=tyme,all=all,fem=fem,male=male)

    return (tyme, all)

if __name__ == "__xmain__":

     yearly()
    
def hold():
    climatology(1969,1979) # will skip 1979
    climatology(1980,1989) # will skip 1982, 1986
    climatology(1990,1999)  
    climatology(2000,2010)  
     

