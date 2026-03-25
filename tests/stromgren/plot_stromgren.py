import numpy as np
import yt

def alpha_b():
    return 2.54e-13 #*(T/1e4)**(-0.8163-0.0208*np.log10(T/1e4))

def rec_time(nH):
    return 1/nH/alpha_b()

def r_strom(n, Q):
    r = (3*Q/4/np.pi/n**2/alpha_b())**(1./3.)
    return r

def r_type(r_s, n, t):
    return r_s*(1-np.exp(-t/rec_time(n)))**(1./3.)*yt.units.cm

def d_type(r_s, c, t):
    r_i = r_s*(1.0+(7./4.)*np.sqrt(4./3.)*c*t/r_s)**(4./7.)
    return r_i*yt.units.cm

def find_snapshots(data_dir, sim_name):
    from glob import glob 
    import os
    pattern = os.path.join(data_dir, f"{sim_name}_hdf5_plt_cnt_*")
    files = sorted(glob(pattern))
    return files

def find_stromgren_radius(filename):
    ds = yt.load(filename)
    ad = ds.all_data()

    sp = ds.sphere(ds.domain_center, ds.domain_right_edge[0])
    rp = yt.create_profile(
        sp,
        ("index", "radius"),
        ("flash", "ihp "),
        units={("index", "radius"): "pc"},
        logs={("index", "radius"): False},
    )

    sim_rs = rp.x.value[np.argmin(abs(rp["flash", "ihp "].value - 0.5))]
    sim_time = ds.current_time.in_units("Myr")

    return sim_time, sim_rs

if __name__=="__main__":

    import matplotlib.pyplot as plt
    # plot analytic solution
    # TODO: user just inputs Q & n, test and automates everything else
    t = np.linspace(0,1,100)*yt.units.Myr
    Q = 1e48 # 1/s
    n = 100 # 1/cm^3
    T = 1e4 # K
    c = 8e5 # cm/s, cs = sqrt(k_B*T/mu/mH), mu=1.3 (hard-coded in torch), T=1e4K

    plt.plot(t, d_type(r_strom(n,Q),c, t.to('s').value).to('pc'), color='r', label='D-type')
    plt.plot(t, r_type(r_strom(n,Q),n, t.to('s').value).to('pc'), color='k', label='R-type')

    for s_d in find_snapshots('data', 'dtype'):
        sim_time, sim_rs = find_stromgren_radius(s_d)
        plt.scatter(sim_time, sim_rs, c='r', marker='x')
    for s_r in find_snapshots('data', 'rtype'):
        sim_time, sim_rs = find_stromgren_radius(s_r)
        plt.scatter(sim_time, sim_rs, c='k', marker='x')

    plt.legend(frameon=False)
    plt.ylabel('r (pc)')
    plt.xlabel('t (Myr)')
    plt.xlim(0,1.0)
    plt.ylim(0,8)
    plt.savefig('stromgren_test.png')
