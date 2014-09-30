#!/usr/bin/env python

'''
Copyright (c) 2014, Kyle Simek
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

'''
Converts the svg file output by the Ellipse and Line Segment Detector (elsd) into
a text tile containing sampled points along the curves.  

Usage: elsd_to_curveset.py <input_svg> <sampling_delta> <output_text_file>

Output format:

    <num curves>
    <curve1 num_pts>
    <curve1 pt1 x> <curve1 pt1 y>
    <curve1 pt2 x> <curve1 pt2 y>
    ...
    <curveN num_pts>
    <curveN pt1 x> <curveN pt1 y>
    ...
'''

from sys import argv
import re
import numpy as np
import math

def lrange(r1, inc, r2):
    n = ((r2-r1)+2*np.spacing(r2-r1))//inc
    return np.append(np.linspace(r1,r1+inc*n,n+1), r2)

def parse_header_(f):
    line = ""
    pattern = re.compile("<svg width=\"([0-9]+)px\" height=\"([0-9]+)px\"")
    m = None
    while True:
        line = f.readline()
        if not line:
            print "Failed to parse header"
            exit(1)

        m = pattern.match(line)
        if m is not None:
            break
    f.readline() # discard second line

    return (m.group(1), m.group(2))

def decimate_line_(L, delta):
    x1, y1, x2, y2 = L

    x0 = np.array([x1, y1])
    x1 = np.array([x2, y2])
    d = x1 - x0
    magnitude = np.linalg.norm(d)
    ts = lrange(0, delta, magnitude)

    tmp = np.array([x0 + d * t/magnitude for t in ts])

    return tmp

def angle_between(u, v):
    theta = math.acos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
    if (u[0] * v[1] - u[1] * v[0]) < 0:
        theta = -theta
    return theta
        
def ellipse_endpoint_to_center(P):
    """
    Convert between SVG-style elliptical path to centered representation
    """
    x1, y1, rx, ry, phi, fa, fs, x2, y2 = P

    phi = math.radians(phi)

    # M is transposed to operate on row vectors
    M = np.array([[math.cos(phi), -math.sin(phi)],[math.sin(phi), math.cos(phi)]])
    v = np.array([(x1 - x2)/2, (y1 - y2)/2])
    x_ = np.dot(v, M)
    rx_sq = rx*rx
    ry_sq = ry*ry
    x1_ = x_[0]
    y1_ = x_[1]
    x1_sq = x1_ * x1_
    y1_sq = y1_ * y1_

    num = rx_sq * ry_sq - rx_sq * y1_sq - ry_sq * x1_sq
    den = rx_sq * y1_sq + ry_sq * x1_sq
    c_ = math.sqrt(num / den) * np.array([rx*y1_/ry, -ry*x1_/rx])

    if fa == fs:
        c_ = -c_

    cx_ = c_[0]
    cy_ = c_[1]

    c = np.dot(c_, M.T) + np.array([(x1+x2)/2, (y1+y2)/2])

    x_hat = np.array([1, 0])
    v1 = np.array([(x1_ - cx_)/rx, (y1_ - cy_)/ry])
    v2 = np.array([(-x1_ - cx_)/rx, (-y1_ - cy_)/ry])
    theta_1 = angle_between(x_hat, v1)
    delta_theta = angle_between(v1, v2)
    
    return (c[0], c[1], rx, ry, phi, theta_1, delta_theta)

def decimate_centered_ellipse_(E, angle_delta):
    cx, cy, rx, ry, phi, theta_1, delta_theta = E

    # M is transposed to operate on row vectors
    M = np.array([[math.cos(phi), math.sin(phi)],[-math.sin(phi), math.cos(phi)]])

    if delta_theta < 0:
        thetas = lrange(theta_1, -angle_delta, theta_1 + delta_theta)
    else:
        thetas = lrange(theta_1, angle_delta, theta_1 + delta_theta)

    out = np.array([ [math.cos(theta), math.sin(theta)] for theta in thetas])

    out[:,0] = rx * out[:,0]
    out[:,1] = ry * out[:,1]

    out = np.dot(out, M)
    out[:,0] += cx
    out[:,1] += cy

    return out
 
def chord_length_parameterize(pts):
    tmp = np.diff(pts, 1, 0)
    return np.cumsum(np.append(0, np.sum(tmp**2, axis=-1))**0.5)

def resample_curve(pts, delta):
    t = chord_length_parameterize(pts)
    assert(len(t) == len(pts))
    t_out = lrange(0, delta, t[-1])
    x_out = np.interp(t_out, t, pts[:,0])
    y_out = np.interp(t_out, t, pts[:,1])



    tmp = np.array([x_out, y_out]).T
    return tmp
    
def decimate_ellipse_path_(P, delta):
    E = ellipse_endpoint_to_center(P)
    cx, cy, rx, ry, phi, theta_1, delta_theta = E

    rad_delta = delta / max(rx, ry)

    pts = decimate_centered_ellipse_(E, rad_delta)
    return resample_curve(pts, delta)

def decimate_ellipse_(E, delta):
    x1, y1, phi, rx, ry = E

    phi = math.radians(phi)

    theta_1 = 0
    delta_theta = 2*math.pi
    E = (x1, y1, rx, ry, phi, theta_1, delta_theta)

    rad_delta = delta / max(rx, ry)
    pts = decimate_centered_ellipse_(E, rad_delta)
    return resample_curve(pts, delta)

def decimate_circle_(C, delta):
    cx, cy, rx, ry = C
    E = (cx, cy, 0, rx, ry)
    return decimate_ellipse_(E, delta)

def parse_geometry_string_(string, delta):
    if string.startswith("<line"):
        pattern = re.compile("<line x1=\"([^\"]+)\" y1=\"([^\"]+)\" x2=\"([^\"]+)\" y2=\"([^\"]+)\"")
        m = pattern.match(string)
        assert(m is not None)
        g = m.groups()
        assert(len(g) == 4)
        L = map(float, g);
        return decimate_line_(L, delta)
    elif string.startswith("<path"):
        pattern = re.compile("<path d=\"M ([^ ,]+),([^ ,]+) A([^ ,]+),([^ ,]+) ([^ ,]+) ([^ ,]+),([^ ,]+) ([^ ,]+),([^ ,]+)\"")
        m = pattern.match(string)
        assert(m is not None)
        g = m.groups()
        assert(len(g) == 9)
        P = map(float, g)
        return decimate_ellipse_path_(P, delta)
    elif string.startswith("<ellipse transform"):
        pattern = re.compile("<ellipse transform=\"translate\(([^ ()\"]+) ([^ ()\"]+)\) rotate\(([^ ()\"]+)\)\" rx=\"([^ ()\"]+)\" ry=\"([^ ()\"]+)\" ")
        m = pattern.match(string)
        assert(m is not None)
        g = m.groups()
        assert(len(g) == 5)
        E = map(float, g)
        return decimate_ellipse_(E, delta)
    elif string.startswith("<ellipse cx"):
        pattern = re.compile("<ellipse cx=\"([^\"]+)\" cy=\"([^\"]+)\" rx=\"([^\"]+)\" ry=\"([^\"]+)\"")

        m = pattern.match(string)
        assert(m is not None)
        g = m.groups()
        assert(len(g) == 4)
        C = map(float, g)
        return decimate_circle_(C, delta)
    else:
        print "Couldn't parse line: \n"  + string
        exit(1)
    
def write_curves_(fname, curves, raw):
    """ 
    raw: set to true to leave out headers (only a long list of points is written, useful for reading and plotting in matlab)"
    """
    f = open(fname, 'w')
    if not raw:
        f.write("%d\n" % len(curves))
    for c in curves:
        if c is None:
            continue
        if not raw:
            f.write("%d\n" % len(c))
        for row in c:
            assert(len(row) == 2)
            f.write("%f %f\n" % (row[0], row[1]))
    f.close()

def main():
    if len(argv) < 4:
        print "usage: elsd_to_curveset.py infile.svg delta out.curves"
    in_fname = argv[1]
    delta = float(argv[2])
    out_fname = argv[3]
    f = open(in_fname)

    [width,height] = parse_header_(f)

    curves = []
    for line in f:
        if line.strip() == "</svg>":
            break
        curves.append(parse_geometry_string_(line, delta))

    raw = False
    write_curves_(out_fname, curves, raw)
    
if __name__ == "__main__":
    main()
