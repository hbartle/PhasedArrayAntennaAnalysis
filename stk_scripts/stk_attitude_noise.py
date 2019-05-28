#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.





import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
#plt.style.use('classic')
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams.update({'font.size': 14})
import matplotlib.gridspec as gridspec
import pandas
import datetime

from util import *


"""
    Apply Attitude noise to a STK Attitude File for a given scenario
"""


def euler_to_quaternion(yaw, pitch, roll):
    qx = np.sin(roll / 2) * np.cos(pitch / 2) * np.cos(yaw / 2) - np.cos(roll / 2) * np.sin(pitch / 2) * np.sin(yaw / 2)
    qy = np.cos(roll / 2) * np.sin(pitch / 2) * np.cos(yaw / 2) + np.sin(roll / 2) * np.cos(pitch / 2) * np.sin(yaw / 2)
    qz = np.cos(roll / 2) * np.cos(pitch / 2) * np.sin(yaw / 2) - np.sin(roll / 2) * np.sin(pitch / 2) * np.cos(yaw / 2)
    qw = np.cos(roll / 2) * np.cos(pitch / 2) * np.cos(yaw / 2) + np.sin(roll / 2) * np.sin(pitch / 2) * np.sin(yaw / 2)

    return np.asarray([qx, qy, qz, qw]).T

def quaternion_to_euler(q):
        x = q[:,0]
        y = q[:,1]
        z = q[:,2]
        w = q[:,3]

        t0 = +2.0 * (w * x + y * z)
        t1 = +1.0 - 2.0 * (x * x + y * y)
        X = np.arctan2(t0, t1)

        t2 = +2.0 * (w * y - z * x)
        #t2[np.where(t2 > 1.0)] = 1.0
        #t2[np.where(t2 < -1.0)] = -1.0
        Y = np.arcsin(t2)

        t3 = +2.0 * (w * z + x * y)
        t4 = +1.0 - 2.0 * (y * y + z * z)
        Z = np.arctan2(t3, t4)

        return np.asarray([X, Y, Z]).T

def quaternion_multiply(q1,q2):
    x1 = q1[:,0]
    y1 = q1[:,1]
    z1 = q1[:,2]
    w1 = q1[:,3]

    x2 = q2[:,0]
    y2 = q2[:,1]
    z2 = q2[:,2]
    w2 = q2[:,3]

    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2

    return np.asarray([x, y, z, w]).T

def quaternion_normalize(q):
    x = q[:,0]
    y = q[:,1]
    z = q[:,2]
    w = q[:,3]

    norm = np.sqrt(np.power(x,2) + np.power(y,2) + np.power(z,2) + np.power(w,2)).reshape(-1,1)
    norm = np.tile(norm,(1,4))
    return q/norm

def fftnoise(f):
    f = np.array(f, dtype='complex')
    Np = (len(f) - 1) // 2
    phases = np.random.rand(Np) * 2 * np.pi
    phases = np.cos(phases) + 1j * np.sin(phases)
    f[1:Np+1] *= phases
    f[-1:-1-Np:-1] = np.conj(f[1:Np+1])
    return np.fft.ifft(f).real

def band_limited_noise(min_freq, max_freq, samples=1024, samplerate=1):
    freqs = np.abs(np.fft.fftfreq(samples, 1/samplerate))
    f = np.zeros(samples)
    idx = np.where(np.logical_and(freqs>=min_freq, freqs<=max_freq))[0]
    f[idx] = 1
    return fftnoise(f)

if __name__ == "__main__":
    output_directory = "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/"

    std = [1,2,3,4,5,6,7,8,9,10]
    std = [30]

    for target_steering_noise_std in std:

        filenames = ["../../stk/Scenario_SAT_SAT_attitude_noise/LeadingSat.a",
               "../../stk/Scenario_SAT_SAT_attitude_noise/TrailingSat.a"]

        target_filenames = ["../../stk/Scenario_SAT_SAT_attitude_noise/LeadingSatNoiseSigma"
                            +str(target_steering_noise_std) + ".a",
                            "../../stk/Scenario_SAT_SAT_attitude_noise/TrailingSatNoiseSigma"
                            +str(target_steering_noise_std) + ".a"]

        for i,filename in enumerate(filenames):
            target_filename = target_filenames[i]

            data = np.genfromtxt(filename, skip_header=25,skip_footer=2)
            timestamp = data[:, 0]
            quaternions = data[:, 1:]

            dt = np.diff(timestamp).mean()
            #print(dt)
            upper_frequency = 0.1
            noise_x = band_limited_noise(0,upper_frequency,samples=timestamp.size,samplerate=1/dt)
            noise_x = target_steering_noise_std*noise_x/noise_x.std()
            noise_y = band_limited_noise(0,upper_frequency,samples=timestamp.size,samplerate=1/dt)
            noise_y = target_steering_noise_std*noise_y/noise_y.std()

            #plt.plot(timestamp,noise_x)
            #plt.plot(timestamp,noise_y)

            noise_z = np.zeros(timestamp.shape)


            noise = euler_to_quaternion(deg2rad(noise_x),deg2rad(noise_y),deg2rad(noise_z))
            noise = noise.reshape((-1,4))

            quaternions = quaternion_normalize(quaternions)
            noise = quaternion_normalize(noise)

            quaternions_noisy = quaternion_multiply(noise,quaternions)
            quaternions_noisy = quaternion_normalize(quaternions_noisy)

            euler = quaternion_to_euler(quaternions)
            euler_noise = quaternion_to_euler(quaternions_noisy)

            # Plot
            gs = gridspec.GridSpec(3, 1)
            fig = plt.figure()
            fig.canvas.set_window_title("Sigma=" +str(target_steering_noise_std))
            ax0 = plt.subplot(gs[0])
            ax1 = plt.subplot(gs[1], sharex=ax0)
            ax2 = plt.subplot(gs[2], sharex=ax0)
            end = 1000

            ax0.plot(timestamp[:end], rad2deg(euler_noise[:end,0]),label="Noise")
            ax0.plot(timestamp[:end], rad2deg(euler[:end,0]), label="No Noise")

            ax0.set_ylabel("Roll [deg]")
            ax0.grid(True)

            ax1.plot(timestamp[:end], rad2deg(euler_noise[:end,1]),label="Noise")
            ax1.plot(timestamp[:end], rad2deg(euler[:end,1]), label="No Noise")

            ax1.set_ylabel("Pitch [deg]")
            ax1.grid(True)

            ax2.plot(timestamp[:end], rad2deg(euler_noise[:end,2]),label="Noise")
            ax2.plot(timestamp[:end], rad2deg(euler[:end,2]), label="No Noise")

            ax2.set_ylabel("Yaw [deg]")
            ax2.set_xlabel("Time [s]")
            ax2.grid(True)

            dic = {"Time": timestamp[:end],
                   "RollNoise": rad2deg(euler_noise[:end,0]),
                   "RollNoNoise": rad2deg(euler[:end,0]),
                   "PitchNoise": rad2deg(euler_noise[:end,1]),
                   "PitchNoNoise": rad2deg(euler[:end,1]),
                   "YawNoise": rad2deg(euler_noise[:end, 2]),
                   "YawNoNoise": rad2deg(euler[:end, 2])}

            df = pandas.DataFrame(data=dic)
            df.to_csv(output_directory + "attitude_noise_sigma="+str(target_steering_noise_std) + ".txt",
                      index=None, sep=",", float_format='%.2f')

            #ax1.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %Y '))
            #ax1.xaxis.set_major_locator(mdates.DayLocator())
            #plt.gcf().autofmt_xdate()
            plt.legend(loc="upper right")

            data_out = np.vstack([timestamp,quaternions_noisy.T]).T

            with open(filename) as sourcefile:
                header = [next(sourcefile) for x in range(25)]

            with open(target_filename,'w') as targetfile:
                targetfile.write(''.join(header))
            with open(target_filename,'ab') as targetfile:
                np.savetxt(targetfile,data_out)
                targetfile.write("\nEND Attitude".encode())

    plt.show()


