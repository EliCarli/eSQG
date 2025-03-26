####################################################################################
# #######################  PAPERS TO BE REFERRED TO #################################
# ##
# ##  [1] Lapeyre, G., P. Klein, 2006a: Dynamics of the upper oceanic layers in terms 
# ##      of surface quasigeostrophic theory. J.P.O. 36, 165-176.
# ##
# ##  [2] Klein, P., Hua B.L., G. Lapeyre, X. Capet, S. LeGentil and H. Sasaki., 2008. 
# ##      Upper Ocean Dynamics from High 3-D Resolution Simulations.J.P.O., 38, 8, 1748???1763. 
# ##
# ##  [3] Klein P., J. Isern-Fontanet, G. Lapeyre, G. Roullet, E. Danioux, 
# ##      B. Chapron, S.Le Gentil and H. Sasaki, 2009. 
# ##      Diagnosis of vertical velocities in the upper ocean from high resolution 
# ##      Sea Surface Height.VOL. 36, L12603, doi:10.1029/2009GL038359. 
# ##
# ##  [4] Carli, E., Siegelman, L., Morrow, R., and Vergara, O., 2024. 
# ##      Surface Quasi Geostrophic Reconstruction of Vertical Velocities and Vertical 
# ##      Heat Fluxes in the Southern Ocean: Perspectives for SWOT. 
# ##      Journal of Geophysical Research: Oceans, 129(9). https://doi.org/10.1029/2024JC021216
# ##
# ##  [5] Carli, E., Tranchant, Y-T., Siegelman, L., Le Guillou, F., Morrow, R., Ballarotta, M., 
# ##      Vergara, O., 2024. Small-scale eddy diagnostics around the Southern Ocean Polar Front 
# ##      with SWOT. ESS Open Archive . January 11, 2025. DOI: 10.22541/essoar.173655546.61867308/v1
# ##
# ################################################################################### 
# ##
# ##      COPYRIGHT: Elisa Carli, Lia Siegelman, Patrice KLEIN, April 2025 
# ##
# ##
# ###################################################################################

import numpy as np



def doubly_per(data):
    """ make doubly periodic by mirror symmetry """
    data = np.hstack((data,np.fliplr(data)))
    data_per = np.vstack((np.flipud(data),data))
    return data_per

def detrend_harmonic(phi):
    '''
    Detrend using the harmonics. Performs worst in our tests
    '''
    M,L = phi.shape
    delta = np.mean(phi[-2,:])-np.mean(phi[1,:])
    a = np.cos(np.pi*np.arange(M)/(M-1))
    a = a[:,None]
    b = np.ones(L)
    b = b[:,None].T
    y = phi - 0.5 * (1-np.dot(a,b)) * delta
    return y

def get_kxky(L,M,dx,dy):
    Lx = L*dx
    Ly = M*dy
    coefx = Lx / (2*np.pi)
    coefy = Ly / (2*np.pi)

    kx = np.hstack((np.arange(L/2),0,np.arange(-L/2+1,0)))/coefx
    ky = np.hstack((np.arange(M/2),0,np.arange(-M/2+1,0)))/coefy

    kkx,kky = np.meshgrid(kx,ky)
    kk = np.sqrt(kkx**2 + kky**2)

    return kx,ky,kkx,kky,kk

def correlation_profile(var1, var2, dim='k'):
    '''
    The two arrays in imput must be xarrays
    '''
    corr = np.zeros(len(var1[dim]))

    for k in np.arange(len(var1[dim])):

        aa_tmp = np.sum(var1.isel(dim=k)*var2.isel(dim=k))
        bb_tmp = np.sum(var1.isel(dim=k)**2)
        cc_tmp = np.sum(var2.isel(dim=k)**2).real
        corr[k] = aa_tmp/ np.sqrt(bb_tmp*cc_tmp)
        
    return corr

def ssh_preprocessing_sqg(ssh):
    
    ssh_doub_periodic = doubly_per(ssh)
    ssh_det1 = detrend_harmonic(ssh_doub_periodic)
    ssh_det1_nomean = ssh_det1 - np.mean(ssh_det1)
    ssh_det2 = detrend_harmonic(ssh_det1_nomean)
    
    return ssh_det2

def ssh_timeseries_preprocessing_sqg(ssh):
    ssh_doub_periodic = doubly_per(ssh)

    for t in range(ssh_doub_periodic.shape[2]):
        ssh_det1 = detrend_harmonic(ssh_doub_periodic[:,:,t])
        ssh_det1_nomean = ssh_det1 - np.mean(ssh_det1)
        ssh_doub_periodic[:,:,t] = detrend_harmonic(ssh_det1_nomean)
    
    return ssh_doub_periodic


## RECONSTRUCTION OF THE VERTICAL VELOCITY FIELD from SSH ##
def sqg_w(ssh, cc, rho_cst, M, L, N, kk, kkx, kky, zz, goN02orho0, gof0):
    
    '''
    M and L are the size of the double periodic field
    At the end of the reconstruction cuts 3 pixels at each border of each direction (x and y) to remove some spectral residuals
    '''
    
    spec_u_ssh=-1j*gof0*kky*np.fft.fft2(ssh)
    spec_v_ssh= 1j*gof0*kkx*np.fft.fft2(ssh)
    spec_rho_ssh=-(kk*np.fft.fft2(ssh))*rho_cst

    u_ssh=np.fft.ifft2(spec_u_ssh).real
    v_ssh=np.fft.ifft2(spec_v_ssh).real
    rho_ssh=np.fft.ifft2(spec_rho_ssh).real

    jac1s=u_ssh*rho_ssh
    jac2s=v_ssh*rho_ssh

    spec_jac1= 1j*kkx*np.fft.fft2(jac1s)
    spec_jac2= 1j*kky*np.fft.fft2(jac2s)
    spec_jacs= spec_jac1+spec_jac2


    kz=0
    for k in range(0,N): 

        if kz == 0:
            w_sqg_matr = np.zeros([len(np.arange(0,N)),int(M/2-6),int(L/2-6)])

        func_exp = np.zeros((M, L))
        func_exp = np.exp(kk*zz[k])

        v_z= np.fft.ifft2(spec_v_ssh*func_exp).real
        u_z= np.fft.ifft2(spec_u_ssh*func_exp).real
        rho_z= np.fft.ifft2(spec_rho_ssh*func_exp).real

        jac1z=u_z*rho_z
        jac2z=v_z*rho_z
        spec_jac1= 1j*kkx*np.fft.fft2(jac1z)
        spec_jac2= 1j*kky*np.fft.fft2(jac2z)
        spec_jacz= spec_jac1+spec_jac2
        jacz=np.fft.ifft2(spec_jacz).real

        spec_jacsz= spec_jacs*func_exp
        jacs=np.fft.ifft2(spec_jacsz).real

        w_sqg=-jacs + jacz
        w_sqg= goN02orho0*w_sqg*cc
        w_sqg=86400*np.flipud(w_sqg[3:int(np.ceil((M/2)-3)), 3:int(np.ceil((L/2)-3))])
        w_sqg_matr[kz,:,:] = w_sqg

        kz = kz+1
        
    return(w_sqg_matr)


## RECONSTRUCTION OF THE RELATIVE VORTICITY FIELD from SSH ##
def sqg_rel_vort(ssh, M, L, N, kk, zz, gof0):
    
    '''
    M and L are the size of the double periodic field
    At the end of the reconstruction cuts 3 pixels at each border of each direction (x and y) to remove some spectral residuals
    '''

    spec_vortssh = - gof0*(kk**2)*np.fft.fft2(ssh)

    spec_vortssh_z=np.zeros((M,L))
    spec_vortssh1=spec_vortssh

    kz = 0
    for k in range(0,N):
        if kz == 0:
            vort_sqg_matr = np.zeros([len(np.arange(0,N)),int(M/2-6),int(L/2-6)]) 

        print(r'k = '+str(k))

        func_exp = np.zeros((M, L))
        func_exp = np.exp(kk*zz[k])
        spec_vortssh_z=spec_vortssh1*func_exp
        xis_sqg=np.flipud(np.fft.ifft2(spec_vortssh_z).real[3:int(np.ceil((M/2)-3)), 3:int(np.ceil((L/2)-3))])
        vort_sqg_matr[k,:,:] = xis_sqg

        kz = kz+1
        
    return vort_sqg_matr


## RECONSTRUCTION OF THE DENSITY FIELD from SSH ##
def sqg_rho(ssh,M,L,N,kk,rho_cst,cc, zz_rmean):
    
    '''
    Function to compute the SQG density
    M and L are the size of the double periodic field
    At the end of the reconstruction cuts 3 pixels at each border of each direction (x and y) to remove some spectral residuals
    '''
        
    xis=np.zeros((M,L))
    spec_rhoz=np.zeros((M,L))

    spec_rhoss=-kk*np.fft.fft2(ssh)
    spec_rhoss1=spec_rhoss*rho_cst

    kz=0
    for k in range(0,N):
        k_coas = k
        if kz == 0:
            rho_sqg_matr = np.zeros([len(np.arange(0,N)),int(M/2-6),int(L/2-6)])

        func_exp = np.zeros((M, L))
        func_exp = np.exp(kk*zz_rmean[k])
        spec_rhoz = spec_rhoss1*func_exp
        rhos_sqg = np.flipud(np.fft.ifft2(spec_rhoz).real[3:int(np.ceil((M/2)-3)), 3:int(np.ceil((L/2)-3))])*cc
        rho_sqg_matr[kz,:,:] = rhos_sqg

        kz = kz+1
        
    return rho_sqg_matr 
