"""Python cover to Gigatraj Trajectory model (written in C++).

  The main class GIGATRAJ has methods to initialize parcels for
gigatraj, running the trajectory model and plotting the
results. Configuration parameters are contained in an YAML file,
typically "gigatraj.yaml". Campaign specific geometry is contained in
a separate yaml file, e.g., "inspyre.yaml".

Arlindo da Silva, August 2026 

(c) 2026 Asticou Earth Systems, LLC

"""

import os
import numpy as np
import yaml
import xarray as xr

from concurrent.futures import ThreadPoolExecutor, as_completed
import shlex
import subprocess

from datetime import datetime, timedelta
from pathlib  import Path

MAP_REGION = dict ( siberia = (22.,-155,25,80),
                    siberia2 = (140.,-100,35,85),
                    north_america = (-120.0,-70.0,22.5,60.0),
                    north_america2 = (-125.0,-70.0,22.5,60.0),
                    north_america3 = (-150.0,-100.0,22,60.0),
                    alaska = (-170.0,-130.0,40,65.0),
                   )

class TRAJError(Exception):
    """
    Defines general exception errors.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class GIGATRAJ(object):

    def __init__(self,config='gigatraj.yaml'):
        """
        Initialize Gigatraj. 
        """

        with open("gigatraj.yaml", encoding="utf-8") as stream:
            self.cf = yaml.safe_load(stream)

        self.Verbose = self.cf['Verbose']

        # Forecast cycle default
        # ----------------------
        if self.cf['Trajectories']['ForecastCycle'] == 'today' :
            self.cf['Trajectories']['ForecastCycle'] = \
                str(datetime.now().date()).replace('-','')+'_00'
        elif self.cf['Trajectories']['ForecastCycle'] == 'yesterday' :
            back = timedelta(days=1)
            self.cf['Trajectories']['ForecastCycle'] = \
                str((datetime.now()-back).date()).replace('-','')+'_00'

        # Map region
        # ----------
        region = self.cf['Plots']['BoundingBox']
        if region in MAP_REGION:
            self.cf['Plots']['Region'] = region
            self.cf['Plots']['BoundingBox'] = MAP_REGION[region]
        else:
            self.cf['Plots']['Region'] = None
            
    def genParcels(self):
        """
        Generate parcels around each fire to start trajectories from.
        """

        cf, F, P, T = self.cf, self.cf['Fires'], self.cf['Parcels'], self.cf['Trajectories']

        # Create list of fires being released
        # -----------------------------------
        fires = []
        for r in cf['Trajectories']['Releases']:
                fires += r['Fires']
 
        # Create parcels for each unique fire
        # -----------------------------------
        for f in list(dict.fromkeys(fires)):

            lon, lat, dz = F[f]['lon'], F[f]['lat'], P['Vdelta_km']
            
            for z in F[f]['altitudes_km']:

                filename = f"{cf['SCRATCH']}/parcel_init.{f}.{z}km.nc" # omitting release time

                zlow, zhigh = z - dz/2, z + dz/2
            
                cmd = f"{cf['PREFIX']}/bin/gt_generate_parcels --random " + \
                      f" --clat {lat} --clon {lon} "     + \
                      f" --zlow {zlow} --zhigh {zhigh} " + \
                      f" --number {P['Number']} " + \
                      f" --radius {P['Radius_km']} --vertical {P['Vertical']} --vunits {P['Vunits']} " + \
                      f" --format {P['Format']} --netcdf {filename}"

                print(cmd+'\n')
                if os.system(cmd):
                    print("*** Error running case ***")
            
    def _prepTrajectories(self):
        """
        Prepares a list of commands to be executed later for generating
        Trajectories given parcels on files produced by genParcels().
        """

        Cmd = []
        
        cf, F, P, T = (self.cf, self.cf['Fires'],
                       self.cf['Parcels'], self.cf['Trajectories'])

        timestep = T['Timestep_min']
        source = T['MetSource']
        catalog = T['Catalog']
        MetSpec = f"ModelRun={T['ForecastCycle']}"

        dirname = f"{cf['OUTPUT']}/{catalog}/data/{T['ForecastCycle']}"
        if os.system(f"mkdir -p {dirname}"):
            raise TRAJError(f"Could not create directory {dirname}")

        # Loop over parcel releases
        # -------------------------
        for R in cf['Trajectories']['Releases']:

            Start, ZeroTime = R['Start'], R['Start']
            Stop = Start + timedelta(days=T['Duration_days'])

            # Create parcels for each unique fire
            # -----------------------------------
            for f in R['Fires']:

                for z in F[f]['altitudes_km']:

                    parcels = f"{cf['SCRATCH']}/parcel_init.{f}.{z}km.nc" # omitting release time 
                    filename = f"{dirname}/parcel_traj.{T['ForecastCycle']}.{f}.{Start.isoformat()[:-3]}_{z}km.nc"
                
###                          f" --metoptions catlog={catalog}.cat;{MetSpec}" + \
###                          f" --metoptions {MetSpec}" + \
                    cmd = f"{cf['PREFIX']}/bin/gtmodel_s01 --verbose" +\
                          f" --begdate {Start.isoformat()} --enddate {Stop.isoformat()}" +\
                          f" --zerodate {ZeroTime.isoformat()}" +\
                          f" --tstep {timestep}" +\
                          f" --source {source}" + \
                          f" --metoptions {MetSpec}" + \
                          f" --vertical {T['VTrace']}" +\
                          f" --parcels {parcels}" +\
                          f" --parcelvertical  {P['Vertical']}" +\
                          f" --input_netcdf " + \
                          f" --frequency {T['outFreq_hour']}" + \
                          f" --netcdf_out {filename}" \
                          f" --format \'{T['OutFormat']}\'"

                    Cmd += [cmd,]

        # All done
        # --------
        return Cmd
                    
    def _listTrajFiles(self, start=None, fire=None):
        """
        Prepares a list of files for all releases.
        """

        myFiles = dict()
        myAlts  = dict()
        
        cf, F, P, T = (self.cf, self.cf['Fires'],
                      self.cf['Parcels'], self.cf['Trajectories'])

        source = T['MetSource']
        catalog = T['Catalog']

        # Loop over parcel releases
        # -------------------------
        for R in cf['Trajectories']['Releases']:

            s = R['Start']

            myFiles[s] = dict()
            myAlts[s] = dict()
            dirname = f"{cf['OUTPUT']}/{catalog}/data/{T['ForecastCycle']}"

            # For each file
            # -------------
            for f in R['Fires']:

                myFiles[s][f] = []
                myAlts[s][f] = []
                
                for z in F[f]['altitudes_km']:

                    filename = f"{dirname}/parcel_traj.{T['ForecastCycle']}.{f}.{s.isoformat()[:-3]}_{z}km.nc"
                    
                    myFiles[s][f] += [filename,]
                    myAlts[s][f]  += [z,]

        # Optionally trim list of files
        # -----------------------------
        if start is not None:
            for s in list(myFiles.keys()):
                if s != start:
                    del myFiles[s] # drop this start
                    del myAlts[s]  # drop this start

        nfires = 0
        if fire is not None:
            for s in list(myFiles.keys()):
                 for f in list(myFiles[s].keys()):
                     if f != fire:
                         #print(f"Dropping <{f}>")
                         del myFiles[s][f] # drop this fire
                         del myAlts[s][f]  # drop this fire
                     else:
                         nfires += 1

            if nfires < 1:
                raise RuntimeError(f"No trajectory file found for fire {fire}")

        return myFiles, myAlts

                    
    def genTrajectories(self):
        """
        Generate Trajectories given parcels on files produced by genParcels().
        """

        # Symlink the met catalog as MetGEOSfp.cat
        # I really don't like this. Apparently we should be
        # able to provide the catalog with --metoptions
        # but I could not get it to work.
        # Using a lock to prevent 2 simultaneous processes
        # to overwrite this.
        # -------------------------------------------------
        lock = 'MetGEOSfp.lock'
        catalog = self.cf['Trajectories']['Catalog']
        if Path(lock).exists():
            raise TRAJError(f"Found file {lock}; if there is no other gigatraj process running, remove it and resubmit")
        else:
            cmd = f"touch MetGEOSfp.lock; /bin/rm -f MetGEOSfp.cat; /bin/ln -s {catalog}.cat MetGEOSfp.cat"
            if os.system(cmd):
                raise TRAJError(f"Cannot execute {cmd}")
                    
        # Generate list of commands to be executed in parallel
        # ----------------------------------------------------
        Cmd = self._prepTrajectories()

        maxConcurrent = self.cf['Trajectories']['MaxConcurrent']
        dry_run = self.cf['Trajectories']['DryRun']
        
        # Execute commands in parallel
        # ----------------------------
        parallel_execute(Cmd, dry_run=dry_run, \
                         Verbose=self.Verbose, maxConcurrent=maxConcurrent)

        # Remove lock when done
        # ---------------------
        cmd = f"/bin/rm -f {lock}"
        if os.system(cmd):
            raise TRAJError(f"Cannot execute {cmd}")

        
    def loadTrajectories(self, start=None, fire=None, merge_fires=False):
        """
        For a given Release *start* datetime, it lazy loads all netcdf files
        associated with this release date. If fire=None it loads all fires,
        otherwise, only files for the specified fire.
        Of course, the trajectories are assumed to have been calculated
        already with method genTrajectories().

        Optionally it can merge all fires into a single merged set of
        Trajectories. When merging, all fires MUST have the same number
        of altitudes.
        """

        # Get a list of files associated with this run
        # --------------------------------------------
        myFiles, myAlts = self._listTrajFiles(start,fire)
        
        # Lazy load the files, annotating with release altitude
        # -----------------------------------------------------
        Trajs = [] # flat list of Trajectories for each start/fire
        for s in myFiles:

            # Load datasets into dictionary
            # -----------------------------
            T_ = dict() 
            for f in myFiles[s]:
                T_[f] = dict()     # T_[f][z]
                for fn, z in zip(myFiles[s][f],myAlts[s][f]):
                    ds = xr.open_dataset(fn,engine='netcdf4')
                    ds.attrs['altitude_km'] = z
                    ds.attrs['fire'] = f
                    ds.attrs['catalog'] = self.cf['Trajectories']['Catalog']
                    T_[f][z] = ds

            # Add to list, merging fires if required.
            # This is a list of list, each element having
            # a list of datasets with all altitudes.
            # --------------------------------------------
            if merge_fires:
                _T = transpose_dict(T_) # _T[z][f]
                trajs = []
                for z in _T:
                    ds = xr.concat(list(_T[z].values()), dim="id",data_vars="minimal")
                    ds["id"] = range(len(ds.id)) # reset ids.
                    ds.attrs['fire'] = 'Merged'
                    trajs += [ds,]
                Trajs += [trajs,]
            else:
                for f in T_:   # T_[f][z]
                    trajs = []
                    for z in T_[f]:
                        trajs += [T_[f][z],] 
                    Trajs += [trajs,]
                
        return Trajs

    def plotTrajectories(self, start=None, fire=None):
        """
        For a given Release *start* datetime, it lazy loads all netcdf files
        associated with this release date. If fire=None it loads all fires,
        otherwise, only files for the specified fire.
        Of course, the trajectories are assumed to have been calculated
        already with method genTrajectories(). 
        """
        from traj_plot import plot_traj

        cf, T = self.cf, self.cf['Trajectories']

        source = T['MetSource']
        catalog = T['Catalog']
        
        dirname = f"{cf['OUTPUT']}/{catalog}/images/{T['ForecastCycle']}"
        region = self.cf['Plots']['Region']
        if region is not None:
            dirname = f"{dirname}/{region}" # different directories for each region
        
        if os.system(f"mkdir -p {dirname}"):
            raise TRAJError(f"Could not create directory {dirname}")

        # Get a list of files associated with this run
        # --------------------------------------------
        myFiles, myAlts = self._listTrajFiles(start,fire)

        # Bounding box
        # ------------
        bbox = cf['Plots']['BoundingBox']
        sats = cf['Plots']['Satellites']
        
        # Loop over and plot
        # ------------------
        #dt = timedelta(hours=self.cf['Plots']['Density_timestep_hours'])
        ValidTime = self.cf['Plots']['Density_valid_time']
        CampaignFile = cf['Plots']['CampaignFile']
        dpi = cf['Plots']['Image_dpi']
        for s in myFiles:
            for f in myFiles[s]:

                # Bundle for this start, fire
                # ---------------------------
                Traj = [] # Will hold multiple release altitudes
                for fn, z in zip(myFiles[s][f],myAlts[s][f]):
                    ds = xr.open_dataset(fn,engine='netcdf4')
                    ds.attrs['altitude_km'] = z
                    ds.attrs['fire'] = f
                    ds.attrs['catalog'] = self.cf['Trajectories']['Catalog']
                    Traj += [ds,]

                StartTime = datetime.fromisoformat(Traj[0].attrs['Trajectory_start'])
                #ValidTime = StartTime + dt
                t0, tv = StartTime.isoformat()[:-3], ValidTime.isoformat()[:-3], 
 
                fig, _ = plot_traj (Traj, ValidTime, CampaignFile,
                                    satellites=sats,
                                    geographic_bounds=bbox)

                img_file = f"{dirname}/parcel_traj.density.{f}.{t0}+{tv}.png"
                if cf['Verbose']:
                    print(f"[] Saving image {img_file}") 

                fig.savefig(img_file, dpi=dpi)
 
        return 
    
    
#-----
def run_cmd(cmd):
    """
    Run a single command as a subprocess.
    """
    command = shlex.split(cmd)
    result = subprocess.run(command,
                            text=True,
                            capture_output=True,
                            check=True)
    return cmd, result.stdout

def parallel_execute(Cmd,dry_run=False,Verbose=False,maxConcurrent=None):
    """
    Run list of commands in Cmd in parallel with ThreadPools.
    """

    # Dry run
    # -------
    if dry_run:
        for cmd in Cmd:
            print(cmd)
        return

    if maxConcurrent is None:
        gt = GIGATRAJ()
        maxConcurrent = gt.cf['Trajectories']['MaxConcurrent']
        Verbose = gt.Verbose
        
    # Serial Run without ThreadPools
    # ------------------------------
    if maxConcurrent == 0:
        for cmd in Cmd:
            print(cmd+'\n')
            if os.system(cmd):
                print(f"*** Error running {cmd} ***")
        return
    
    # Normal Concurrent Run
    # ---------------------
    with ThreadPoolExecutor(max_workers=maxConcurrent) as executor:

        futures = [executor.submit(run_cmd, cmd) for cmd in Cmd]

        for future in as_completed(futures):

            try:
                cmd, stdOut = future.result()
                if Verbose:
                    print('[] '+cmd,'\n',stdOut)

            except subprocess.CalledProcessError as error:
                print(f"Command failed: {error}")

            except Exception as error:
                print(f"Unexpected error: {error}")

    print("Parallel execution completed.")

def transpose_dict(d):
    """
    If "d" is a directory addressed as d[a][b] it returns
    a dictionary addressed as d_T[b][a]. It is assumed that
    the inner keys "b" are uniform.
    """
    # Extract unique inner keys (b)
    inner_keys = {b for inner in d.values() for b in inner}

    return {b: {a: d[a][b] for a in d if b in d[a]} for b in inner_keys}
#--------------------------------------------------------------------------------

if __name__ == "__main__":

    gt = GIGATRAJ()

    #myFiles = gt.loadReleases()
