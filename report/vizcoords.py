# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 12:40:23 2021

@author: slulla
"""


def adjust_coordinates(exon_coords, var_coords, intron_size=10):
    
    '''
    UNIT TEST 1 - intron at 0
    >>> ec1 = [(100,200), (250,400), (700,900)]
    >>> vc1 = [150, 300, 350, 750, 850]
    >>> adjust_coordinates(ec1, vc1)
    ([(10, 110), (120, 270), (280, 480)], [60, 170, 220, 330, 430])
    
    UNIT TEST 2 - exon at 0
    >>> ec2 = [(0,200), (250,400), (700,900)]
    >>> vc2 = [150, 350, 750, 850]
    >>> adjust_coordinates(ec2, vc2)
    ([(0, 200), (210, 360), (370, 570)], [150, 310, 420, 520])
    
    UNIT TEST 3 - intron at 0, intron_size=20
    >>> ec3 = [(100,200), (250,400), (700,900)]
    >>> vc3 = [150, 300, 350, 750, 850]
    >>> adjust_coordinates(ec3, vc3, intron_size=20)
    ([(20, 120), (140, 290), (310, 510)], [70, 190, 240, 360, 460])
    
    UNIT TEST 4 - exon at 0, intron_size=20
    >>> ec4 = [(0,200), (250,400), (700,900)]
    >>> vc4 = [150, 350, 750, 850]
    >>> adjust_coordinates(ec4, vc4, intron_size=20)
    ([(0, 200), (220, 370), (390, 590)], [150, 320, 440, 540])
    
    UNIT TEST 5 - 1 intron length < intron_size
    >>> ec5 = [(0,200), (250,400), (410,610)]
    >>> vc5 = [150, 350, 460, 510]
    >>> adjust_coordinates(ec5, vc5, intron_size=20)
    ([(0, 200), (220, 370), (380, 580)], [150, 320, 430, 480])

    UNIT TEST 5 - 1 intron length < intron_size
    >>> ec5 = [(0,200), (250,400), (410,610)]
    >>> vc5 = [150, 350, 460, 510]
    >>> adjust_coordinates(ec5, vc5, intron_size=20)
    ([(0, 200), (220, 370), (380, 580)], [150, 320, 430, 480])

    UNIT TEST 7 - variation coordinates inside introns
    >>> ec7 = [(0,200), (250,400), (700,900)]
    >>> vc7 = [150, 225, 350, 480, 750, 850, 1000]
    >>> adjust_coordinates(ec7, vc7)
    ([(0, 200), (210, 360), (370, 570)], [150, 310, 420, 520])
    
    UNIT TEST 8 - backwards coordinates
    >>> ec8 = [(900,700), (400,250), (200,100)]
    >>> vc8 = [850, 750, 350, 300, 150]
    >>> adjust_coordinates(ec8, vc8)
    ([(10, 110), (120, 270), (280, 480)], [60, 170, 220, 330, 430])

    UNIT TEST 9 - variation coordinate == start or end
    >>> ec9 = [(100,200), (300,350), (400,600)]
    >>> vc9 = [150, 200, 400, 500]
    >>> adjust_coordinates(ec9, vc9)
    ([(10, 110), (120, 170), (180, 380)], [60, 110, 180, 280])

    UNIT TEST 10 - 1 intron, 1 exon, empty var_coords
    >>> ec10 = [(100,200)]
    >>> vc10 = []
    >>> adjust_coordinates(ec10, vc10)
    ([(10, 110)], [])

    UNIT TEST 11 - empty input
    >>> ec11 = []
    >>> vc11 = []
    >>> adjust_coordinates(ec11, vc11)
    ([], [])

    UNIT TEST 12 - overlapping exon coordinates
    >>> ec12 = [(0,200), (100,300), (400,500)]
    >>> vc12 = [100, 200, 300]
    >>> adjust_coordinates(ec12, vc12)
    Traceback (most recent call last):
        ...
    ValueError: Exon coordinates overlap
    
    

    
    '''
    
    
    #param exon_coords: [(int,int)]; backwards exon coords
    #return flipped_exon_coords: [(int,int)]; forwards exon coords
    def flip_exon_coords(exon_coords):
        flipped_exon_coords = []
        for (start,end) in exon_coords[::-1]: #flip exon_coords without permanently changing it
            flipped_exon_coords.append((end,start))
        return flipped_exon_coords
    
    
        
    #if empty exon_coords, return empty exon coords and original var coords
    if len(exon_coords) == 0: return [],var_coords
    
    #if exon coordinates are backwards, flip them
    if exon_coords[0][0] > exon_coords[0][1]: exon_coords = flip_exon_coords(exon_coords)
    
    compact_exon_coords = []
    compact_var_coords = []
    end_last_compact = 0 #end of last exon in compacted coords
    end_last_original = 0 #end of last exon in original coords
    
    for (start,end) in exon_coords:
        
        if end_last_original > start: raise ValueError('Exon coordinates overlap')
        
        if start==0:
            compact_exon_coords.append((start,end))
            end_last_compact = end
            end_last_original = end
            compact_var_coords = list(filter(lambda x: x < end,var_coords))
            continue
        
        intron_len = start - end_last_original #length of intron between this exon and previous
        new_start = end_last_compact+intron_size if intron_len > intron_size else end_last_compact+intron_len
        new_end = new_start + end - start ###1 based or 0 based coordinate system?###
        compact_exon_coords.append((new_start, new_end))
        
        #find all variation coords corresponding to this exon
        for c in list(filter(lambda x: x>=start and x<=end, var_coords)):
            compact_var_coords.append(new_start + c - start)
        
        end_last_compact = new_end
        end_last_original = end
        
#     print(compact_exon_coords)
#     print(sorted(compact_var_coords))
    
    #sort variation coordinates ascending
    return compact_exon_coords,sorted(compact_var_coords)

#param exons: [(int,int)]; list of exon coordinates
#return flat_exons: [(int,int)]; list of exon coords where overlapping coords are condensed into longest exon possible
def flatten_exons(exons):
    
    '''
    UNIT TEST 1 - no overlaps
    >>> e1 = [(150, 200), (0,100), (300, 450), (500, 600)]
    >>> flatten_exons(e1)
    [(0, 100), (150, 200), (300, 450), (500, 600)]
    
    UNIT TEST 2 - all overlaps
    >>> e2 = [(50, 100), (75, 150), (0, 80), (120, 130)]
    >>> flatten_exons(e2)
    [(0, 150)]

    UNIT TEST 3 - start with overlaps
    >>> e3 = [(10, 50), (40, 70), (60, 80), (100, 150), (160, 200), (300, 350)]
    >>> flatten_exons(e3)
    [(10, 80), (100, 150), (160, 200), (300, 350)]

    UNIT TEST 4 - end with overlaps
    >>> e4 = [(50, 100), (150, 200), (250, 350), (375, 390), (325, 400), (270, 320), (400, 450), (380, 450), (450, 500)]
    >>> flatten_exons(e4)
    [(50, 100), (150, 200), (250, 500)]

    UNIT TEST 5 - overlaps only in middle
    >>> e5 = [(200, 250), (275, 300), (350, 400), (310, 400), (360, 380), (400, 410), (450, 500)]
    >>> flatten_exons(e5)
    [(200, 250), (275, 300), (310, 410), (450, 500)]
    
    '''
    
    def largest_range(ls):
        min_start = min([start for (start,end) in ls])
        max_end = max([end for (start,end) in ls])
        return (min_start, max_end)

    exons = sorted(exons, key=lambda x:x[0]) #sort by starting coord
    exons_flat = [] #list to return
    overlaps = [exons[0]] #list of overlapping exons
    max_end = exons[0][1] #end of overlapping region

    for i in range(1, len(exons)): 
        if exons[i][0] <= max_end: #if start of this exon is before end of overlapping region:
            overlaps.append(exons[i]) #add this exon to overlapping region
            if exons[i][1] > max_end: max_end = exons[i][1]
        else:
            if len(overlaps) == 1: exons_flat.append(overlaps[0])
            if len(overlaps) > 1: exons_flat.append(largest_range(overlaps))
            overlaps = [exons[i]]
            max_end = exons[i][1]

    exons_flat.append(largest_range(overlaps)) #have to do this one more time at the end
    exons_flat = sorted(exons_flat, key=lambda x:x[0]) #sort by starting coord
    
    return exons_flat

def var_coord_in_exon_coords(var_coord, exon_coords):
    for (start,end) in exon_coords:
        if var_coord >= start and var_coord <= end: return True
    return False



if __name__ == "__main__":
    import doctest
    doctest.testmod()