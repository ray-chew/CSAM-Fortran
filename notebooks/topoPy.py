import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import scipy.signal as sig
from sklearn.linear_model import LinearRegression
import matplotlib as mpl

##########################################################################################

def readnc(datfile,var):
    #---- topography data file
    df = nc.Dataset(datfile)
    print(df.variables)
    
    #---- variable and coordinates
    lon, lat = df.variables['lon'][:], df.variables['lat'][:]
    z = df.variables[var['name']][:]
    df.close()
    
    return lon, lat, z

##########################################################################################

def fetch_data(datfile,var,lon_centre,lat_centre,lon_width,lat_width):
    
    #---- variable and the coordinates
    lon, lat, z = readnc(datfile, var)
    
    #---- get number of records,nlat,nlon
    nrecords = np.shape(z)[0]; nlon = np.shape(lon)[1]; nlat = np.shape(lat)[1]
    
    #---- process each record to get the (lon,lat) for each topographic observation in 1D
    lon_res=[]; lat_res=[]; z_res=[]   # resulting lon,lat,z  
    
    for n in range(nrecords):
        print('n = ',n)
        lon_,lat_ = np.meshgrid(lon[n][:],lat[n][:])
#         print(lon_, lat_)
        lon_= lon_.ravel() 
        lat_ = lat_.ravel() 
        z_ = np.flipud(z[n][:]).ravel() 
#         print(lon_.shape, lat_.shape, z_.shape)
        idx = np.nonzero((np.abs(lon_ - lon_centre)<= lon_width/2) & 
                         (np.abs(lat_ - lat_centre)<= lat_width/2))[0]
        print(idx)
        if len(idx)!=0:
            lon_dummy,lat_dummy,z_dummy = lon_[idx],lat_[idx],z_[idx]
            lon_res.extend(lon_dummy.tolist())
            lat_res.extend(lat_dummy.tolist())
            z_res.extend(z_dummy.tolist())
    
    lon_res = np.array(lon_res)
    lat_res = np.array(lat_res)
    z_res = np.array(z_res)
        
    del lat, lon, z
        
    #---- processing of the lat,lon,topo to get the regular 2D grid for topography
    lon_uniq, lat_uniq = np.unique(lon_res), np.unique(lat_res) # get unique values of lon,lat
    nla = len(lat_uniq); nlo = len(lon_uniq)
    
#     print("lat_res shape = ", lat_res.shape)
#     print("lon_res shape = ", lon_res.shape)
#     print("z_res shape = ", z_res.shape)
#     print("nla = ", nla)
#     print("nlo = ", nlo)

    #---- building 2D topography field
    lat_lon_topo = np.vstack((lat_res,lon_res,z_res)).T
    lat_lon_topo = lat_lon_topo[lat_lon_topo[:,0].argsort()]  # sorted according to latitude
    lon_sort_id = [lat_lon_topo[n*nlo:(n+1)*nlo,1].argsort()+nlo*n for n in range(nla)] 
    lon_sort_id = np.array(lon_sort_id).reshape(-1)
    lat_lon_topo = lat_lon_topo[lon_sort_id]  # sorted according to longitude for each len(lon_u)
    topo_2D = np.reshape(lat_lon_topo[:,2],(nla,nlo))
    del lat_lon_topo, lon_sort_id
        
    print('Data fetched...')
    return lon_uniq, lat_uniq, topo_2D

##############################################################################################

def fft_2D(topo,res_x,res_y): # topo must be 2D
    topospec = np.fft.fft2(topo) # 2D FFT
    topospec = np.fft.fftshift(topospec) # centre for the frequencies
    
    #---- get the wavenumbers
    k_x = np.fft.fftfreq(np.shape(topo)[1], d=res_x)
    k_y = np.fft.fftfreq(np.shape(topo)[0], d=res_y)
    k_x, k_y = np.fft.fftshift(k_x), np.fft.fftshift(k_y)
    
    return k_x,k_y,topospec

###############################################################################################

#---- 2D Fourier spectrum fitting ----#
#---- Approach 1 ----#

def fitFourierSpec(data,Ni,Nj,I,J,nhar_i,nhar_j):
    # data, I, J must be in 1D
    # Ni, Nj: total number of points in i and j indices
    # nhar_i, nhar_j: number of harmonics in i and j index
    
    m_i, n_j = np.arange(-nhar_i+1,nhar_i), np.arange(-nhar_j+1,nhar_j)
    dof = len(m_i)*len(n_j)

    # basis matrix
    basis = []
    for k in range(len(data)):
        coeff = [np.exp(1j*2*np.pi*(mm*I[k]/Ni + nn*J[k]/Nj)) for mm in m_i for nn in n_j]
        basis.append(coeff)

    basis = np.array(basis)
    #print(basis.shape)

    # Solve: data_k = \sum basis_kl*a_l + \epsilon
    # obtain a_l using least square minimization of the error \epsilon

    h_tilda_l = np.zeros((dof,), dtype='complex')
    E_tilda_lm = np.zeros((dof,dof), dtype='complex')

    for l in range(dof):
        h_tilda_l[l] = np.sum(data*basis[:,l])
        E_tilda_lm[l,:] = np.sum(basis*np.expand_dims(basis[:,l], axis=1), axis=0)

    # now invert E_tilda_lm to get the coefficients
    a_m = np.linalg.inv(E_tilda_lm).dot(h_tilda_l)

    # regular FFT considers normalization by total number of datapoints
    # so multiply the Fourier coefficients by N and return
    return a_m*len(data)

################################################################################################

#---- 2D Fourier spectrum fitting ----#
#---- Approach 2 ----#

def fitFourierSeries(data,Ni,Nj,I,J,nhar_i,nhar_j):
    # data, I, J must be in 1D
    # Ni, Nj: total number of points in i and j indices
    # nhar_i, nhar_j: number of harmonics in i and j index
    
    # number of harmonics for x,y
    m_i, n_j = range(nhar_i), range(nhar_j)

    # basis matrix
    basis_cos = []
    basis_sin = []
    for k in range(len(data)):
        coeff_cos = [np.cos(2*np.pi*(mm*I[k]/Ni + nn*J[k]/Nj)) for mm in m_i for nn in n_j]
        coeff_sin = [np.sin(2*np.pi*(mm*I[k]/Ni + nn*J[k]/Nj)) for mm in m_i for nn in n_j if mm!=0 or nn!=0]
        basis_cos.append(coeff_cos)
        basis_sin.append(coeff_sin)

    coeff = np.hstack([basis_cos,basis_sin])
    tot_coeff = coeff.shape[1]
    print('Total coefficients:',tot_coeff)

    # Solve: data_k = \sum basis_kl*a_l + \epsilon
    # obtain a_l using least square minimization of the error \epsilon
    # alternative

    h_tilda_l = np.zeros((tot_coeff,))
    E_tilda_lm = np.zeros((tot_coeff,tot_coeff))

    for l in range(tot_coeff):
        h_tilda_l[l] = np.sum(data*coeff[:,l])
        E_tilda_lm[l,:] = np.sum(coeff*np.expand_dims(coeff[:,l], axis=1), axis=0)

    # now invert E_tilda_lm to get the coefficients
    a_m = np.linalg.inv(E_tilda_lm).dot(h_tilda_l)

    # regular FFT considers normalization by total number of datapoints N=100
    # so multiply the Fourier coefficients by N here
    a_m = a_m*len(data)

    mid = (a_m.size+1)//2
    fourier_coeff = (a_m[1:mid] + 1j*a_m[mid:])/2    # half the amplitudes are the Fourier coefficients
    fourier_coeff = np.insert(fourier_coeff,0,a_m[0])
    fourier_coeff = fourier_coeff.reshape((nhar_i,nhar_j))
    
    # reconstruct the dataset
    data_recons = coeff.dot(a_m)/len(data)
    
    return fourier_coeff, data_recons

##################################################################################################

#---- return the vector given two points 
def vector(x1,y1,x2,y2):
    return [x2-x1, y2-y1]

#---- determine if a point is inside the triangle
def pointInTriangle(vx,vy,x,y):
    # Inputs:
    # vx: x-coordinate of the vertices
    # vy: y-coordinate of the vertices
    # x,y: coordinates of the target point
    #
    # Output:
    # Boolean: point is inside or outside
    
    x1, x2, x3 = vx
    y1, y2, y3 = vy
    e1 = vector(x1,y1,x2,y2) # edge 1
    e2 = vector(x2,y2,x3,y3) # edge 2
    e3 = vector(x3,y3,x1,y1) # edge 3
    
    p2e1 = vector(x,y,x1,y1) # point to edge 1
    p2e2 = vector(x,y,x2,y2) # point to edge 2
    p2e3 = vector(x,y,x3,y3) # point to edge 3
    
    c1 = np.cross(e1,p2e1)  # cross product 1
    c2 = np.cross(e2,p2e2)  # cross product 2
    c3 = np.cross(e3,p2e3)  # cross product 3
    
    return np.sign(c1) == np.sign(c2) == np.sign(c3)

####################################################################################################
