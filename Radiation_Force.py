# -*- coding: utf-8 -*-
"""
Time: 2020/Jul/06

@author: Benchi Zhao
Python version: 3.7
Functionality: This module calculate the force caused by radiation pressure.
Dependence: Need import 'NewNonABCD', all parameters of the experiment is saved in that file 'input'.
"""

import sys
import numpy as np
import input
import NewNonABCD as NonABCD
import copy


all_forwards = []
all_backwards = []
def intensity(position):
    I = []
    for i in range(len(position)):
        intens = NonABCD.gaussian(position[i], input.sigma)
        I.append(intens)
    normalise_I = I / sum(I)
    return normalise_I

def bundle():
    list_of_rays = np.linspace(-input.width, input.width, input.no_of_rays)
    for i in range(input.no_of_rays):
        NonABCD.ray(0, list_of_rays[i], 0, intensity(list_of_rays)[i])
        NonABCD.lens_trace(input.lens_f, input.lens_pos, input.len_thickness)
        NonABCD.free_propagate(10)
        NonABCD.propagate(input.droplet_pos[0] - input.lens_pos - 10)
        NonABCD.circle(input.radius)
        NonABCD.free_propagate(15)
        NonABCD.propagate(200)
        a = copy.deepcopy(NonABCD.forward)
        b = copy.deepcopy(NonABCD.backward)
        NonABCD.forward.clear()
        NonABCD.backward.clear()
        all_forwards.append(a)
        all_backwards.append(b)

def useful_data_for():
    '''
    Clean the data, some rays in all_forwards will not interact with the droplet, which are useless.
    This function will keep those rays interact with the droplet and delete those rays have no interaction with the droplet.
    '''
    useful_for = copy.deepcopy(all_forwards)
    for i in range(len(useful_for)-1,-1,-1): # iterate from right to left
        if len(useful_for[i]) != 12:
            useful_for.pop(i)
    return useful_for

def useful_data_back():
    '''
    Clean the data, some rays in all_backwards will not interact with the droplet, which are useless.
    This function will keep those rays interact with the droplet and delete those rays have no interaction with the droplet.
    '''
    useful_back = copy.deepcopy(all_backwards)
    for i in range(len(useful_back)-1,-1,-1): # iterate from right to left
        if len(useful_back[i]) != 12:
            useful_back.pop(i)
    return useful_back

def F_s():
    '''
    Calculate the scattering force.
    :return total_forcr: list
        The total_foce contains two element, the first one is the force along x-aixs, and the second one is the force along y-axis.
    '''
    data = useful_data_for()
    force = []
    total_force = np.array([0, 0])
    for i in range(len(data)):
        x = data[i][6][0]
        y = data[i][6][1]
        tilt = np.arcsin((y - input.droplet_pos[1]) / input.radius)
        theta1 = data[i][6][2] - tilt  # Incident angle
        theta2 = NonABCD.snell(theta1)  # Refracted angle
        P = input.power * data[i][6][3]
        F = input.medium_n * P / input.c * (1+NonABCD.R(2*theta1)*np.cos(2*theta1)- (NonABCD.T(theta1) ** 2 * (np.cos(2 * theta1 - 2 * theta2) + NonABCD.R(theta1) * np.cos(2 * theta1))) / (1 + NonABCD.R(theta1) ** 2 + 2 * NonABCD.R(theta1) * np.cos(2 * theta2)))
        Fx = abs(F) * np.cos(theta1)
        Fy = abs(F) * np.sin(theta1)
        # print(Fx,Fy,F)
        force.append([Fx, Fy])
        total_force = total_force + np.array([Fx, Fy])
        # print(total_force)
        # print(' ')
    return total_force

def F_g():
    '''
        Calculate the gradient force.
        :return total_forcr: list
            The total_foce contains two element, the first one is the force along x-aixs, and the second one is the force along y-axis.
        '''
    data = useful_data_for()
    force = []
    total_force = np.array([0, 0])
    for i in range(len(data)):
        x = data[i][6][0]
        y = data[i][6][1]
        tilt = np.arcsin((y-input.droplet_pos[1])/input.radius)
        theta1 = data[i][6][2]-tilt # Incident angle
        theta2 = NonABCD.snell(theta1) # Refracted angle
        P = input.power * data[i][6][3]
        F = input.medium_n*P/input.c * (NonABCD.R(theta1)*np.sin(2*theta1)-(NonABCD.T(theta1)**2*(np.sin(2*theta1-2*theta2)+NonABCD.R(theta1)*np.sin(2*theta1)))/(1+NonABCD.R(theta1)**2+2*NonABCD.R(theta1)*np.cos(2*theta2)))
        Fx = abs(F) * np.cos(np.pi-tilt)
        Fy = abs(F) * np.sin(np.pi-tilt)
        # print(Fx, Fy, F)
        force.append([Fx,Fy])
        total_force = total_force + np.array([Fx, Fy])
        # print(total_force)
    return total_force

if __name__ == '__main__':
    bundle()
    useful_data_for()
    useful_data_back()
    F_g()
    F_s()
    print('F_g=',F_g())
    print('F_s=',F_s())
    print('total force=',np.array(F_s())+np.array(F_g()))