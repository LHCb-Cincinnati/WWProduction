""" Stuff to add:
    - Make CheckRadiativeDecay Faster.
"""
# Imports
import pdb
from numpy import sqrt

def DeltaR(particle1, particle2):
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

def DeltaRMatchBool(reco_particle, mc_particle, max_deltar_val=0.1):
    delta_r = DeltaR(reco_particle, mc_particle)
    if (delta_r > max_deltar_val):
        return(False)
    else:
        return(True)

def FindDecayProduct(truth_particles, source_pid, num_min_products, target_particles_list):
    source_particles_dict = {particle.pt():particle for particle in truth_particles
                    if (particle.particleID().pid()==source_pid
                        and not CheckRadiativeDecay(particle)
                        and CheckOriginVertex(particle))}
    highest_pT_source = sorted(source_particles_dict.keys(), reverse=True)[0]
    source_particle = source_particles_dict[highest_pT_source]
    source_decay_products_dict = {particle_ref.pt():particle_ref.target() for particle_ref
                        in source_particle.endVertices()[0].products()
                        if CheckMotherDaughterCharge(particle_ref)
                        and (particle_ref.particleID().pid() in target_particles_list)}
    highest_pT_target = source_decay_products_dict[max(source_decay_products_dict.keys())]
    return(highest_pT_target)

def CheckRadiativeDecay(particle):
    particle_pid = particle.particleID().pid()
    decay_products = particle.endVertices()[0].products()
    decay_pid = [particle.particleID().pid() for particle in decay_products]
    if particle_pid in decay_pid:
        return(True)
    else:
        return(False)

def CheckOriginVertex(particle):
    mother = particle.mother()
    if particle.originVertex().type()==1:
        return(True)
    elif not mother:
        return(False)
    else:
        return(CheckOriginVertex(mother))

def CheckMotherDaughterCharge(particle, same_polarity=True):
    mother = particle.mother()
    if not mother:
        return(False)
    mother_charge = mother.particleID().threeCharge()
    daughter_charge = particle.particleID().threeCharge()
    polarity_ratio = daughter_charge / mother_charge
    if polarity_ratio>0:
        return(True)
    else:
        return(False)
