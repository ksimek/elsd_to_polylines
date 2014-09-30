elsd_to_polylines
=================

Simple program to parse the SVG output of the ["Ellipse and Line Segment Detector"][1] and outputs piecewise linear curves (polylines).

Note: This is not an SVG parser; it assumes the specific format of the current (1.0) version of elsd.

Usage: 
---------

  elsd_to_curveset.py <input_svg> <sampling_delta> <output_text_file>

Output format:
-------------

    <num curves>
    <curve1 num_pts>
    <curve1 pt1 x> <curve1 pt1 y>
    <curve1 pt2 x> <curve1 pt2 y>
    ...
    <curveN num_pts>
    <curveN pt1 x> <curveN pt1 y>
    ...

[1]: http://ubee.enseeiht.fr/vision/ELSD/ "Ellipse and Line Segment Detector"
