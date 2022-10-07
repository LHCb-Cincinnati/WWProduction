
def DeltaR(particle1, particle2):
    from math import sqrt

    p1_vec = particle1.momentum()
    p2_vec = particle2.momentum()
    delta_eta = p1_vec.eta() - p2_vec.eta()
    delta_phi = p1_vec.phi() - p2_vec.phi()
    delta_r = sqrt(delta_eta**2 + delta_phi**2)
    return(delta_r)

def FindMCParticles(mc_particles_collection, pid):
    pid_mc_particles = [particle
                        for particle in mc_particles_collection 
                        if abs(particle.particleID().pid()) == pid]
    return(pid_mc_particles)

def DeltaRMatching(reco_particle, mc_particles_collection, max_deltar_val=0.1):
    mc_particles_list = FindMCParticles(mc_particles_collection, 13)
    deltar_list = [DeltaR(mc_particle, reco_particle)
                   for mc_particle in mc_particles_list]
    min_deltar_val = min(deltar_list)
    if (min_deltar_val > max_deltar_val):
        return(False)
    else:
        min_deltar_index = deltar_list.index(min(deltar_list))
        return(mc_particles_list[min_deltar_index])
