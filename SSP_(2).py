import os
from datetime import datetime, timezone
from numpy import pi, cos, sin, log, exp, arccos
import numpy as np
import matplotlib.pyplot as plt
from mycode.position import*
import math

class SSP:
    """A Class for handling Sound Speed Profile data"""

    def __init__(self):

        # The data attributes
        self.obs_time = None
        self.log_time = None
        self.obs_latitude = None
        self.obs_longitude = None
        self.vessel_latitude = None
        self.vessel_longitude = None
        self.obs_sample = list()
        self.obs_depth = list()
        self.obs_ss = list()
        self.proc_ss = np.array([])
        self.twtt_layer=np.array([])
        
        #from LabB
        self.vessel_speed = None
        self.bot_depth = None
        self.d = np.array([])
        self.c = np.array([])
        self.g = np.array([])
       
        
        self.metadata = dict()
        self.metadata["units"] = "rad"
        self.metadata["count"] = None
        self.metadata["geoid"] = None
        self.metadata["ellipsoid"] = None
        self.metadata["chart_datum"] = None
        self.metadata["time_basis"] = "UTC"

    # The I/O methods:
    
    def read_mvp_file(self,fullpath):
        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening sound speed profile data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)
            
        # Open, read and close the file
        svp_file = open(fullpath)
        svp_content = svp_file.read()
        svp_file.close    
        
        # Tokenize the contents
        svp_lines = svp_content.splitlines()

        # Create a Position object       
        #pos = Position()
        
        
        # extract the header lines
        h_lines = 0
        for line in svp_lines:
            #print(line)
            if (len(line)) ==0: # stop reading file after for space
                break
            h_lines += 1
        
#         print(h_lines) # number of lines taken by header
                
        for line in svp_lines:
            # parse time
            if "GPS Time" in line:
                obs = line.split()
                #break
                
                self.obs_time = ParseNMEA0183_ZDA(obs[2]) 
                             
            #parse position
            #for line in svp_lines:
            if "GPS Position" in line:
                obs = line.split()
                #break
                
                gga = ParseNMEA0183_GGA(obs[2]) 
                #print(gga)
                self.obs_latitude = gga[1]
                self.obs_longitude = gga[2]
              
            # Parse the Depth
            #for line in svp_lines:
            if "Bottom Depth" in line:
                obs = line.split()
                self.bot_depth=float(obs[2])
        
            # Parse the Vessel Speed
            #for line in svp_lines:
            if "Ship Speed" in line:
                obs = line.split()
                #print(line)
                self.vessel_speed=float(obs[2])/1.852
            
                
        # extract the units for the observations
        rec_type = (svp_lines[h_lines+1].split(',')) # +1 is to read the line after the empty space
        #print(rec_type)
        
        index_depth = rec_type.index('Dpth(m)')
        index_ss = rec_type.index('SV(m/s)')
        
        #print(index_depth)
            
        #parsing the records
        for line in svp_lines[h_lines+2: ]:
            obs_1 = line.split(',')
            #print(obs[3])
            self.obs_depth.append(float(obs_1[index_depth]))
            self.obs_ss.append(float(obs_1[index_ss]))
        
        #sort tthe records (smallest to largest    
        temp = sorted(zip(self.obs_depth, self.obs_ss), key=lambda x: x[0])
        self.obs_depth, self.obs_ss = map(list, zip(*temp))   
            
        # Remove any duplicate depths with associated sound speeds
        d_p=self.obs_depth[0]
        index = 0
        unwanted = []
        for d in self.obs_depth[1:]:
            index += 1
            if d == d_p:
                unwanted.append(index)
            d_p = d

        for e in sorted( unwanted, reverse = True):
            del self.obs_depth[e]
            del self.obs_ss[e]   
            
        # B.0.2.9 Creating numpy arrays
        self.d = np.array(self.obs_depth)
        self.c = np.array(self.obs_ss)
        #self.g = np.array(self.obs_gradient)
                
        #B.0.2.10 Extending the profiles to the Surface
        if self.d[0] > 0:
            self.d = np.insert(self.d,0,0)
            self.c = np.insert(self.c,0,self.c[0])         
            
        #B.0.2.11 Calculating the Sound Speed Gradients
        self.g = (self.c[1:] - self.c[0:-1])/(self.d[1:] - self.d[0:-1])
        
        dd = (self.c[1:] - self.c[0:-1])
        #print("dd %s:" % dd)
        
        #Expanding the profiles to Full Ocean Depth
        if self.d[-1] < 12000:
            x = 0.017*(12000 - self.d[-1]) + self.c[-1]
            self.c = np.append(self.c,[x])
            self.d = np.append(self.d,[12000])
            self.g = np.append(self.g,([x]-self.c[-2])/(12000 - self.d[-2]))
                 
        #B.0.2.13 Replacing Zero Gradients
        self.g[self.g == 0] = 0.00001
        
        #print(self.g.shape, self.d.shape, self.c.shape)
        
        
    def ray_trace(self, d_start, th_start, twtt):
        pass
    
    def determine_depth(self, d_start, th_start, ss_start, twtt):
        #d_start: depth start
        # B1.3.0.0 Initialization
        depth=0
        rad_dist=0
        layer_s=0
        layer_e=0

        # B.1.3.0.1 Determine the start layer
        layer_s = sum(d_start >= self.d) -1
        #print(layer_s) # current layer
        
        # B.1.3.0.2 Determine the Ray Constant
        c_ray = cos(th_start)/ss_start #cos theeta d/sCstart
       # print(c_ray)

        # B.1.3.0.3 Calculate Ray Path Properties for all layer
        r_curve = -(1/(c_ray*self.g[0:])) #ray curve
        #print(r_curve.shape)
        

        #  B1.3.0.0 Return values as tuple

        th = np.arccos(self.c[0:]*c_ray) #depression angle at top of layer for all
        #print(th.shape)
        
#        # print(d_angle)
#         rad_dist = R*(sin(math.acos(self.c[layer_s+1]*c_ray))-sin(d_angle)) # Radius of path curvature in layer i
        dx = r_curve*(sin(th[1:]) - sin(th[0:-1]))
        #print(dx[0:100]) 
        #print(dx.shape)
        rad_dist = np.sum(dx[0:-1])
        
            
        #print(dx.shape)
        
       
        #vert_dist = self.d[layer_s + 1] - self.d[layer_s] #vertical distance traversed in layer
         
        # vertical distance
        dz = self.d[1:] - self.d[0:-1]
        
        #harmonic mean
        hm = (1/self.g) - log(self.c[1:]/self.c[0:-1])
        print(hm)

        # One way travel time
        #owtt = H*((1 + cos(d_angle))/(1+cos(self.c[layer_s+1])
        
        upper = 1 + cos(th[0:-1])
        lower = 1+ cos(th[1: ])
        owtt = hm * upper/lower #hm harmonic mean
        twtt = 2*owtt
        
        print(twtt)
        #twtt = twtt/owtt
         
        #thickness for each layer
        delta_d = (self.d[1:] - self.d[0:-1])
        #print(sum(delta_d[28:29]))
        
        return depth, rad_dist, layer_s, layer_e;  
    
    def read_jhc_file(self, fullpath):

        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening sound speed profile data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        motion_file = open(fullpath)
        motion_content = motion_file.read()
        motion_file.close

        # Tokenize the contents
        motion_lines = motion_content.splitlines()
        self.obs_time = datetime.fromtimestamp(
            float(motion_lines[1].split()[0]), timezone.utc)
        self.log_time = datetime.fromtimestamp(
            float(motion_lines[2].split()[0]), timezone.utc)
        self.obs_latitude = float(motion_lines[3].split()[0])
        self.obs_longitude = float(motion_lines[3].split()[1])
        self.vessel_latitude = float(motion_lines[4].split()[0])
        self.vessel_longitude = float(motion_lines[4].split()[1])
        self.metadata["count"] = int(motion_lines[5].split()[0])

        count = 0  # initialize the counter for the number of rows read

        for motion_line in motion_lines[16:]:
            observations = motion_line.split()  # Tokenize the stringS
            self.obs_sample.append(float(observations[0]))
            self.obs_depth.append(float(observations[1]))
            self.obs_ss.append(float(observations[2]))
            count += 1

        if self.metadata["count"] != count:
            raise RuntimeError('Nr of Samples read ('+str(count) +
                               ') does not match metadata count (' +
                               str(self.metadata["count"])+')')

        # Process the data - in the jhc data files this is already a one-way profile,
        # this just for illustration
        self.proc_ss = np.zeros((count, 3))

        # Sort the data samples by depth
        sorted_ss = sorted(zip(self.obs_depth, self.obs_ss))

        layer = 0
        for d, ss in sorted_ss:
            self.proc_ss[[layer], [0]] = d
            self.proc_ss[[layer], [1]] = ss
            layer += 1

        # Identify all the depths for which there are multiple observations
        mask = np.full((count, 1), True)
        mask[1:, [0]] = np.diff(self.proc_ss[:, [0]], axis=0) != 0

        # Remove the duplicates - You really should get statistical representations here
        # but to keep this short just remove the duplicates
        self.proc_ss = self.proc_ss[mask[:, 0], ...]

        # Determine the gradients - Note the indexing: the gradient of the first layer 
        # is contained at the same index as the data for the TOP of the layer.
        self.proc_ss[0:-1, [2]] = np.diff(self.proc_ss[:, [1]],
                                          axis=0)/np.diff(self.proc_ss[:, [0]], axis=0)

        # Estimate gradient for last layer assuming that the temperature and salinity remain the same
        # gradient solely a function of pressure (depth)
        self.proc_ss[-1, [2]] = 0.017

        # Extend to 12000 m if necesarry - this is to get around some manufcturers requirements
        if self.obs_depth[-1] < 12000:
            ss = self.proc_ss[-1:, [1]] + self.proc_ss[-1:, [2]] \
             * (12000-self.proc_ss[-1:, [0]])
            self.proc_ss = np.vstack((self.proc_ss, [12000, ss, 0.017]))

        # Extend to 0 m if necesarry - assume well mixed
        if self.obs_depth[0] > 0:
            self.proc_ss = np.vstack(
                ([0, self.proc_ss[0, [1]], 0.], self.proc_ss))
            
        # Step 5 Create a look-up array of twtts for each full layer
        # Allows for great gain in efficiency (do not have to calculate for each ping)
        self.twtt_layer = np.zeros((count, 1))
        
        for layer in range(0,self.metadata["count"]-1):
            if self.proc_ss[layer, [2]] == 0:
                self.twtt_layer[layer] = 2 * \
                    (self.proc_ss[layer+1, [0]] - self.proc_ss[layer, [0]])/ \
                     self.proc_ss[layer, [1]]
            else:
                self.twtt_layer[layer] = 2 / self.proc_ss[layer, [2]] * \
                 log(self.proc_ss[layer+1, [1]]/self.proc_ss[layer, [1]])

def ParseNMEA0183_ZDA(dt_str):
    obs = dt_str.split(',')
    time = datetime( 
            int( obs[4]), int( obs[3]), int( obs[2]), 
            int( obs[1][0:2]), int( obs[1][2:4]), int( obs[1][4:6]),int(obs[1][7:])*10000)
    return time

def ParseNMEA0183_GGA(dt_str):

    # Get the GGA string and tokenize it
    gga_data = dt_str.split(',')

    # verify that we have a GGA string
    if not gga_data[0][-3:] == "GGA":
        raise RuntimeError(
                'ParseNMEA0183_GGA: argument `dt_str` must be a GGA message')

    # Determine the time of day from both the header and the GGA string

    gga_timedelta=timedelta(hours=int(gga_data[1][0:2]), \
                             minutes = int(gga_data[1][2:4]), \
                             seconds = int(gga_data[1][4:6]))

    # Parse the latitude
    if gga_data[3].lower() == "n":
        lat = float(gga_data[2][0:2])+float(gga_data[2][2:])/60.
    else:
        lat = -float(gga_data[2][0:2])-float(gga_data[2][2:])/60.             

    # Parse the longitude
    if gga_data[5].lower == "w":
        lon = float(gga_data[4][0:3])+float(gga_data[4][3:])/60.
    else:
        lon = -float(gga_data[4][0:3])-float(gga_data[4][3:])/60.               

    # Parse the GNSS Quality indicator
    q = int(gga_data[6])

    # Parse the number of GNSS satellites used for the solution
    n_sats = int(gga_data[7])

    # Parse the HDOP Quality indicator
    hdop = float(gga_data[8])

    # Parse the orthometric height 
    height = float(gga_data[9])

    # Generate an error if the units of the orthometric height is not meters

    if gga_data[10].lower() != "m":
        raise RuntimeError('Orthomeric height units are not meters!')  

    # Parse the geoid ellipsoid separation
    separation = float(gga_data[9])

    if gga_data[12].lower() != "m":
        raise RuntimeError('Orthomeric height units are not meters!') 

    # If there is more data then parse it
    corr_age = None
    corr_station = None
    if not gga_data[13] == "":
        corr_age = float(gga_data[13])
        corr_station = float(gga_data[14][0:-3])

    return gga_timedelta, lat, lon, q, n_sats, hdop, height, separation, corr_age, corr_station    
    

   