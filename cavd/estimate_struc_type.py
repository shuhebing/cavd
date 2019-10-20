# encoding: utf-8
'''
estimate_struc_type -- Judging the type of atomic packing in a structure

@author:     he bing

@copyright:  2017 School of Computer Engineering and Science, ShangHai University. All rights reserved.


@license:    license

@contact:    bhe@i.shu.edu.cn
@deffield    2017/08/03 Updated
'''

import sys
import os
import time
import xlwt
import math
import pymatgen as pmg
import pymatgen.analysis.structure_matcher as pmgmach
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from pymatgen.core.sites import PeriodicSite
from click.core import batch

__all__ = []
__version__ = 0.1
__date__ = '2017-08-03'
__updated__ = '2018-06-27'

#batch read .cif filename
def batch_read_filename(path,filetype):
    filenames=[]
    for i in os.listdir(path):
        if filetype in i:
            #filenames.append(i.replace(filetype,''))
            filenames.append(i)
    return filenames

#write the result to excel file
def write_to_excel(results):
    #create a workbook,encoding:utf-8
    workbook=xlwt.Workbook(encoding='utf-8')
    #create a sheet
    sheet=workbook.add_sheet('sheet 1')  
    #create style
    style=xlwt.XFStyle()
    #create and set font
    font=xlwt.Font()
    font.name='Times New Roman'
    #add font to style
    style.font=font
    #create alignment and set to center
    alignment = xlwt.Alignment()
    alignment.horz=xlwt.Alignment.HORZ_CENTER
    #add alignment to style
    style.alignment=alignment    
    #create another style style1
    style1 = xlwt.XFStyle()
    #create and set font1
    font1 = xlwt.Font()
    font1.name = 'Times New Roman'
    #Set font colour）
    #font1.colour_index = 3 
    font1.bold = True
    style1.font = font1
    style1.alignment = alignment
    
    #set column name
    columnNames=['filename','bcc','normalized rms displacement and maximum distance between paired sites(bcc)',
                 'fcc','normalized rms displacement and maximum distance between paired sites(fcc)',
                 'hcp','normalized rms displacement and maximum distance between paired sites(hcp)']
    #get column number
    columns=len(columnNames)
    #set column width and name 
    for i in range(columns):
        sheet.col(i).width=6000
        sheet.write(0,i,columnNames[i],style1)
        
    #set row number
    rows=len(results) 
    #insert data
    for i in range(0,rows):
        for j in range(columns):
            sheet.write(i+1,j,results[i].split('\t')[j],style)
    #save excel
    excelTime=time.strftime("%Y%m%d-%H%M%S")
    workbook.save('structure_matcher_result-'+excelTime+'.xls')

#construct template structure bcc, fcc, hcp
def CreateTemplateStructure(strutype,size,specie):
    if strutype=='bcc' :
        lattice=[[size,0,0],[0,size,0],[0,0,size]]
        species=[specie,specie]
        coords=[[0,0,0],[1.0/2.0,1.0/2.0,1.0/2.0]]
        return pmg.Structure(lattice,species,coords)
    elif strutype=='fcc' :
        lattice=[[size,0,0],[0,size,0],[0,0,size]]
        species=[specie]*4
        coords=[[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]]
        return pmg.Structure(lattice,species,coords)
    elif strutype=='hcp' :
        lattice=[[size,0,0],[0,size,0],[0,0,size*math.sqrt(6)*(2/3)]]
        species=[specie]*2
        coords=[[0,0,0],[2.0/3.0,1.0/3.0,1.0/2.0]]
        return pmg.Structure(lattice,species,coords)
    
class cstructure(pmg.Structure): 
    def save_species(self, species):
        """
        Remove all occurrences of several species from a structure.
        
        Args:
            species: Sequence of species to remove, e.g., ["Li", "Na"].
        """
        new_sites = []
        #species = [Specie.from_string(s) for s in species]
    
        for site in self._sites:
            new_sp_occu = {sp: amt for sp, amt in site.species_and_occu.items()
                            if sp.as_dict()['element'] in species}
            if len(new_sp_occu) > 0:
                new_sites.append(PeriodicSite(
                    new_sp_occu, site.frac_coords, self._lattice,
                    properties=site.properties))
        self._sites = new_sites

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg


def struc_matcher_cmd(argv=None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by hebing on %s.
  Copyright 2017 School of Computer Engineering and Science, ShangHai University. All rights reserved.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-e", "--elements", nargs=1,dest="elements",action="store",help=" 堆积原子的元素种类")
        parser.add_argument("-l", "--ltol", dest="ltol",action="store",type=float,default=0.2,help=" 分数长度容差[缺省: %(default)s]")
        parser.add_argument("-s", "--stol", dest="stol",action="store",type=float,default=1,help=" rms长度容差[缺省: %(default)s]")
        parser.add_argument("-a", "--angle_tol", dest="angle_tol",action="store",type=float,default=5,help=" 角度容差（单位：度）[缺省: %(default)s]")
        parser.add_argument("-infile", action="store",help=" 结构文件名")
        parser.add_argument("-batch", action="store",type=bool,default=False,help=" 批量读入cifs文件夹下的结构文件")

        # Process arguments
        args = parser.parse_args()
        lentol=args.ltol
        atol=args.angle_tol
        stol=args.stol
        elements=args.elements
        batch=args.batch
         
        if not batch:
            infile=args.infile
            print(args.infile)
            stru=cstructure.from_file(infile)
            stru.save_species(elements)
            V=stru.volume
            n=stru.num_sites
            stol=0.3*(V/n)**(1/3)

            bccstruct=CreateTemplateStructure('bcc',3.0,stru.species[0])
            fccstruct=CreateTemplateStructure('fcc',3.0,stru.species[0])
            hcpstruct=CreateTemplateStructure('hcp',3.0,stru.species[0])
            
            pcomp=pmgmach.StructureMatcher(ltol=lentol,stol=stol, angle_tol=atol,primitive_cell=False, attempt_supercell=True)
            print('Structure like bcc is {}.'.format(pcomp.fit(stru,bccstruct)),'RMS displacement between two structures and max maximum distance between paired sites is {}'.format(pcomp.get_rms_dist(stru,bccstruct)))
            print('Structure like fcc is {}.'.format(pcomp.fit(stru,fccstruct)),'RMS displacement between two structures and max maximum distance between paired sites is {}'.format(pcomp.get_rms_dist(stru,fccstruct)))
            print('Structure like hcp is {}.'.format(pcomp.fit(stru,hcpstruct)),'RMS displacement between two structures and max maximum distance between paired sites is {}'.format(pcomp.get_rms_dist(stru,hcpstruct)))        
        else:
            path="./cifs/"
            filetype=".cif"
            results=[]
            filenames=batch_read_filename(path, filetype)
            print(filenames)
            for file in filenames:
                infile=path+file
                stru=cstructure.from_file(infile)
                stru.save_species(elements)
                V=stru.volume
                n=stru.num_sites
                stol=0.3*(V/n)**(1/3)
                
                bccstruct=CreateTemplateStructure('bcc',3.0,stru.species[0])
                fccstruct=CreateTemplateStructure('fcc',3.0,stru.species[0])
                hcpstruct=CreateTemplateStructure('hcp',3.0,stru.species[0])
                
                pcomp=pmgmach.StructureMatcher(ltol=lentol,stol=stol, angle_tol=atol,primitive_cell=False, attempt_supercell=True)
                results.append(infile+'\t' \
                               +str(pcomp.fit(stru,bccstruct))+'\t'+str(pcomp.get_rms_dist(stru,bccstruct))+'\t' \
                               +str(pcomp.fit(stru,fccstruct))+'\t'+str(pcomp.get_rms_dist(stru,fccstruct))+'\t' \
                               +str(pcomp.fit(stru,hcpstruct))+'\t'+str(pcomp.get_rms_dist(stru,hcpstruct)))
            #print(results)
            write_to_excel(results)
            print("Matcher completed!")
            
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

def struc_matcher(element, infile, batch=False, ltol=0.2, stol=1, atol=5): # IGNORE:C0111
    try:
        if not batch:
            stru=cstructure.from_file(infile)
            stru.save_species(element)
            V=stru.volume
            n=stru.num_sites
            stol=0.3*(V/n)**(1/3)

            bccstruct=CreateTemplateStructure('bcc',3.0,stru.species[0])
            fccstruct=CreateTemplateStructure('fcc',3.0,stru.species[0])
            hcpstruct=CreateTemplateStructure('hcp',3.0,stru.species[0])
            
            pcomp=pmgmach.StructureMatcher(ltol=ltol,stol=stol, angle_tol=atol,primitive_cell=False, attempt_supercell=True)
            print('Structure like bcc is {}.'.format(pcomp.fit(stru,bccstruct)),'RMS displacement between two structures and max maximum distance between paired sites is {}'.format(pcomp.get_rms_dist(stru,bccstruct)))
            print('Structure like fcc is {}.'.format(pcomp.fit(stru,fccstruct)),'RMS displacement between two structures and max maximum distance between paired sites is {}'.format(pcomp.get_rms_dist(stru,fccstruct)))
            print('Structure like hcp is {}.'.format(pcomp.fit(stru,hcpstruct)),'RMS displacement between two structures and max maximum distance between paired sites is {}'.format(pcomp.get_rms_dist(stru,hcpstruct)))        
        else:
            path=infile
            filetype=".cif"
            results=[]
            filenames=batch_read_filename(path, filetype)
            print(filenames)
            for file in filenames:
                infile=path+file
                stru=cstructure.from_file(infile)
                stru.save_species(element)
                V=stru.volume
                n=stru.num_sites
                stol=0.3*(V/n)**(1/3)
                
                bccstruct=CreateTemplateStructure('bcc',3.0,stru.species[0])
                fccstruct=CreateTemplateStructure('fcc',3.0,stru.species[0])
                hcpstruct=CreateTemplateStructure('hcp',3.0,stru.species[0])
                
                pcomp=pmgmach.StructureMatcher(ltol=ltol,stol=stol, angle_tol=atol,primitive_cell=False, attempt_supercell=True)
                results.append(infile+'\t' \
                               +str(pcomp.fit(stru,bccstruct))+'\t'+str(pcomp.get_rms_dist(stru,bccstruct))+'\t' \
                               +str(pcomp.fit(stru,fccstruct))+'\t'+str(pcomp.get_rms_dist(stru,fccstruct))+'\t' \
                               +str(pcomp.fit(stru,hcpstruct))+'\t'+str(pcomp.get_rms_dist(stru,hcpstruct)))
            #print(results)
            write_to_excel(results)
            print("Matcher completed!")
            
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0