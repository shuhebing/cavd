
'''
Created on 
@author: hebing
'''
import math
import os

from monty.serialization import loadfn
import vtk

import numpy as np
import pymatgen.core.periodic_table  as pmgpp
from vtk.vtkCommonCorePython import vtkStringArray
from vtk.vtkCommonDataModelPython import vtkCellArray
from scipy.spatial import ConvexHull
from numpy import copy
from Structure.Structure import Site,Structure
module_dir = os.path.dirname(os.path.abspath(__file__))
EL_COLORS = loadfn(os.path.join(module_dir, "ElementColorSchemes.yaml"))
class ViewAtomSite(object):
    def __init__(self,parentstruct,datasite=None):

        self._sphere=None
        self._mapper=vtk.vtkOpenGLSphereMapper()
        self._actor=vtk.vtkActor()
        self._ColorMap=vtk.vtkUnsignedCharArray()
        self._ColorMapData=vtk.vtkAppendPolyData()
        self._CartPos=np.array([0,0,0])
        self._colors={}  
        self._atomRadius=0.0
        self._Resolution=16
        self._parentstruct=parentstruct
        self._Visibility=True
        if datasite is not None:
            self.GetDataModelSite(datasite)        
    def GetDataModelSite(self,site=None):
        self.SetSiteId(site.GetSiteId())
        self.SetSiteLabel(site.GetSiteLabel())
        self.SetElementsOccupy(site.GetElementsOccupy())
        self.SetPosition(site.GetPosition())
    def GetColors(self):
        return self._colors
    def SetColors(self,Colors={}):
        for key,value in Colors.items():
            self._colors[key]=value
        self.UpdateDisplay()
    def SetResolution(self,resolution):
        self._Resolution=resolution
        self.__UpdateDisplay()
    def GetResolution(self):
        return self._Resolution
    def SetRadius(self,radius):
        self._atomRadius=radius
        self.__UpdateDisplay()
    def GetRadius(self):
        return self._atomRadius
    def GetActor(self):
        return self._actor
    def GetSiteLabel(self):
        return self._sitelabel
    def SetSiteLabel(self,sitelabel):
        self._sitelabel=sitelabel
    def GetSiteId(self):
        return self._siteId
    def SetSiteId(self,siteid):
        self._siteId=siteid
    def GetElementsOccupy(self):
        return self._Elements
    def SetElementsOccupy(self,eleoccupys):
        self._Elements=eleoccupys
        self._atomRadius=0.0
        occupy=0.0       
        vacexist=False    
        self._colors['Vac']=EL_COLORS["VESTA"]['Vac']
        for (key,value) in self._Elements.items():
            self._colors[key]=EL_COLORS["VESTA"][key]
            if key!='Vac':
                occupy=occupy+value;
                self._atomRadius=self._atomRadius+0.3*pmgpp.Element(key).average_ionic_radius
            else:
                vacexist=True
        self._atomRadius=self._atomRadius/len(self._Elements)
        if occupy<1.0 :
            if vacexist:
                self._Elements['Vac']=self._Elements['Vac']+1-occupy
            else:
                self._Elements['Vac']=1.0-occupy
        self._CreateVis()
    def GetRender(self):
        return self._render
    def GetPosition(self):
        return self._position
    def SetPosition(self,pos=np.array([0,0,0])):
        self._CartPos=self._parentstruct.FracPosToCartPos(pos)
        self.UpdateDisplay()
    def _CreatePartColorSphereData(self):
        
        self._ColorMap.SetNumberOfComponents(3)
        ntups=int(self._Resolution*(self._Resolution/2-2)*2)
        self._ColorMap.SetNumberOfTuples(ntups)
        resphi=self._Resolution/2-1
        res=self._Resolution
        occupy=[]
        colors=[]
        for (key,value) in self._Elements.items():
            occupy.append(value)
            colors.append(self._colors[key])
        ocun=[0,0,0]
        totaloccupy=sum(occupy)
        for (i,ocu) in enumerate(occupy):
            ocun[i]=int(round(ocu*res/totaloccupy))
        if sum(ocun)>res:
            ocun[-1]=ocun[-1]-sum(ocun)+res
        nstart=0
        nend=0
        
        self._ColorMap.Initialize()
        for (i,ocu) in enumerate(colors): 
            nend=nstart+ocun[i]
            nlist= range(nstart,nend)
            for j in nlist:
                self._ColorMap.InsertTuple3(j,colors[i][0],colors[i][1],colors[i][2])   
                self._ColorMap.InsertTuple3(j+res,colors[i][0],colors[i][1],colors[i][2])
                for k in range(int(resphi*2-4)):
                    self._ColorMap.InsertTuple3(int(k+j*(resphi*2-4)+2*res),colors[i][0],colors[i][1],colors[i][2])
            nstart=nend 
    def _CreateDisplaySphere(self):
        self._sphere=vtk.vtkSphereSource()
        #self._sphere.LatLongTessellationOn()
        self._sphere.SetPhiResolution(int(self._Resolution/2))
        self._sphere.SetThetaResolution(self._Resolution)
        self._sphere.SetCenter(self._CartPos)
        self._sphere.SetRadius(self._atomRadius)
        self._sphere.Update()
        self._sphere.GetOutput().GetCellData().SetScalars(self._ColorMap) 
      #  self._ColorMapData.RemoveAllInputs()
        self._ColorMapData.AddInputData(self._sphere.GetOutput())
    def SetVisibility(self,boolvalue):
        self._actor.SetVisibility(boolvalue)
        self._Visibility=boolvalue 
    def GetVisibility(self):
        return self._Visibility
    def SetOpacity(self,opac):
        self._actor.GetProperty().SetOpacity(opac)
    def GetOpacity(self):
        return self._actor.GetProperty().GetOpacity()
  

    def _CreateVis(self):
        self._mapper.SetInputConnection(self._ColorMapData.GetOutputPort());
        self._actor.SetMapper(self._mapper)
        self._parentstruct.GetRender().AddActor(self._actor)
#     def __CreateLabel(self):
#         pts=vtk.vtkPoints()
#         labels=vtk.vtkStringArray()
#         verts=vtk.vtkCellArray()
#         labels.SetName('Labels')
#         pts.InsertNextPoint(self._position)
#         labels.InsertNextValue('test')
#         verts.InsertCellPoint(0)
#         ldm = vtk.vtkPolyDataMapper()
#         pd=vtk.vtkPolyData()
#         pd.SetPoints(pts)
#         pd.SetVerts(verts)
#         pd.GetPointData().AddArray(labels)
# #         hie=vtk.vtkPointSetToLabelHierarchy()
# #         hie.SetInputData(pd)
# #         hie.SetMaximumDepth(15)
# #         hie.SetLabelArrayName('Labels')
# #         hie.SetTargetLabelCount(100)
#   #      hie.SetTextProperty(vtkTextProperty)
#         #strategy = vtk.vtkFreeTypeLabelRenderStrategy()
# #         ids = vtk.vtkIdFilter()
# #         ids.SetInputConnection(pd.GetOutputPort())
# #         ids.PointIdsOn()
# #         ids.CellIdsOn()
# #         ids.FieldDataOn()
#         ldm.SetInputData(pd)
#         #ldm.SetRenderStrategy(strategy)
# #         ldm.SetLabelModeToLabelFieldData()
#     #    ldm.SetFieldDataName('labels')
#         self.__labelactor.SetMapper(ldm)
#         self._render.AddActor(self.__labelactor)
        
    def AddToRender(self,ren=None):
        ren.AddActor(self._actor) 
    def UpdateDisplay(self):
        self._CreatePartColorSphereData() 
        self._CreateDisplaySphere()
    def GetColor(self):
        for (key,color) in self._colors.items():
            if key != 'Vac' :
                return color
class ViewAtomPeriodSite(ViewAtomSite):
    def __init__(self,parent,site=None):
        self._MultiCartpos=None;
        self._spheres=[]
        self._actors=[]
        super(ViewAtomPeriodSite,self).__init__(parent,site)

    def SetPosition(self,pos=np.array([0,0,0])):
        self._position=pos
        self._CartPos=self._parentstruct.FracPosToCartPos(pos)
        self.UpdateDisplay()
    def _CaluPeriodDisplayPos(self):
        expandcell=self._parentstruct.GetExpandCell()
        self._PeriodDisplayPos=[]
        endpos=(expandcell[1]+1.0)
        startpos=expandcell[0]*1.0
        pos=copy(self._position)
        pos=pos+startpos
        while pos[0]>=startpos[0] and pos[0]<=endpos[0]+0.0001:
            pos[1]=self._position[1]+startpos[1]
            while pos[1]>=startpos[1] and pos[1]<=endpos[1]+0.0001:
                pos[2]=self._position[2]+startpos[2]
                while pos[2]>=startpos[2] and pos[2]<=endpos[2]+0.0001:
                    self._PeriodDisplayPos.append(copy(pos))
                    pos[2]+=1.0
                pos[1]+=1.0
            pos[0]+=1.0
    def _CreateDisplaySphere1(self):
#        self._ColorMapData.RemoveAllInputs()
        self._CaluPeriodDisplayPos()
        for pos in self._PeriodDisplayPos:
            sphere=vtk.vtkSphereSource()
            sphere.SetPhiResolution(int(self._Resolution/2))
            sphere.SetThetaResolution(self._Resolution)
            cartpos=self._parentstruct.FracPosToCartPos(pos)
            sphere.SetCenter(cartpos[0],cartpos[1],cartpos[2])
            sphere.SetRadius(self._atomRadius)
            sphere.Update()   
            sphere.GetOutput().GetCellData().SetScalars(self._ColorMap)
            self._ColorMapData.AddInputData(sphere.GetOutput())         
            self._spheres.append(sphere) 
         
    def _CreateVis(self):
        pass   
    def _CreateDisplaySphere(self):
        #self._ColorMapData.RemoveAllInputs()
        self._CaluPeriodDisplayPos()
        for act in self._actors:
            self._parentstruct.GetRender().RemoveActor(act)
        del self._actors[:]
        self._polydatamaps=[]
        self._AppendPolydata=[]
        #self._actors=[]
        for pos in self._PeriodDisplayPos:
            sphere=vtk.vtkSphereSource()
            ator=vtk.vtkFollower()
            atext = vtk.vtkVectorText()
            atext.SetText(self.GetSiteLabel())
            apd=vtk.vtkAppendPolyData()
            self._AppendPolydata.append(apd)
            self._actors.append(ator)
            apd.AddInputData(sphere.GetOutput())
            apd.AddInputData(atext.GetOutput())
            polydatamap=vtk.vtkPolyDataMapper()
            polydatamap.SetInputConnection(apd.GetOutputPort())
            
            self._polydatamaps.append(polydatamap)
            sphere.SetPhiResolution(int(self._Resolution/2))
            sphere.SetThetaResolution(self._Resolution)
            cartpos=self._parentstruct.FracPosToCartPos(pos)
            #sphere.SetCenter(cartpos[0],cartpos[1],cartpos[2])
            ator.SetMapper(polydatamap)
            ator.SetPosition(cartpos)
            sphere.SetRadius(self._atomRadius)
            sphere.Update()   
            sphere.GetOutput().GetCellData().SetScalars(self._ColorMap)
            #self._ColorMapData.AddInputData(sphere.GetOutput())         
            self._spheres.append(sphere)
            self._parentstruct.GetRender().AddActor(ator)
            ator.SetCamera(self._parentstruct.GetRender().GetActiveCamera())
#             ator.GetProperty().SetAmbient(0.1)
#             ator.GetProperty().SetDiffuse(0.6) 
#             ator.GetProperty().SetSpecular(0.7) 
#             ator.GetProperty().SetSpecularPower(3)
            self._actors.append(ator)
    def SetVisibility(self,boolvalue):
        for act in self._actors:
            act.SetVisibility(boolvalue) 
        self._Visibility=boolvalue 
    def SetOpacity(self,opac):
        for act in self._actors:
            act.GetProperty().SetOpacity(opac)
    def GetOpacity(self):
        return self._actors[0].GetProperty().GetOpacity()
        
        
class ViewBond(object):
    def __init__(self,parentstructure,atomstart=np.array([[0,0,0],[0,1,0]]),atomend=np.array([[10,10,10],[0.2,1,1]])):
        self._parent=parentstructure
        self._StartAtomPos=atomstart[0]
        self._DirectVector=atomend[0]-atomstart[0]
        self.__Radius=0.08
        self._Resolution=6
        self.__DisplayType=0
        self.__DisplayColor0=atomstart[1]
        self.__DisplayColor1=atomend[1]
        self.__DisplayOpacity=1.0
        self._Points = vtk.vtkPoints()
        self._line0=vtk.vtkCellArray()
        self._line1=vtk.vtkCellArray()
        self.__pd0= vtk.vtkPolyData()
        self.__pd0.SetPoints(self._Points)
        self.__pd0.SetLines(self._line0)
        self.__pd1= vtk.vtkPolyData()
        self.__pd1.SetPoints(self._Points)
        self.__pd1.SetLines(self._line1)
        
        self.__tube0 = vtk.vtkTubeFilter() 
        self.__tube0.SetInputData(self.__pd0)
        self.__tube0.SetRadius(self.__Radius)
        self.__tube0.SetNumberOfSides(self._Resolution)
        
        self.__tube1 = vtk.vtkTubeFilter() 
        self.__tube1.SetInputData(self.__pd1)
        self.__tube1.SetRadius(self.__Radius)
        self.__tube1.SetNumberOfSides(self._Resolution)
        
        self.__mapper0=vtk.vtkPolyDataMapper()
        self.__mapper1=vtk.vtkPolyDataMapper()
        self.__mapper0.SetInputConnection(self.__tube0.GetOutputPort())
        self.__mapper1.SetInputConnection(self.__tube1.GetOutputPort())
        self.__actor0=vtk.vtkActor()
        self.__actor0.SetMapper(self.__mapper0)
        self.__actor1=vtk.vtkActor()
        self.__actor1.SetMapper(self.__mapper1)
        
    def SetEndPoint(self,startpos,endpos):
        self._StartAtomPos=startpos
        self._DirectVector=endpos-startpos
        midposition=self._StartAtomPos+self._DirectVector*0.5
        cartstartpos=self._parent.FracPosToCartPos(self._StartAtomPos)
        cartmidposition=self._parent.FracPosToCartPos(midposition)
        cartendpos=self._parent.FracPosToCartPos(endpos)
        
        self._Points.InsertPoint(0,cartstartpos)
        self._Points.InsertPoint(1,cartmidposition)
        self._Points.InsertPoint(2,cartendpos)
        
        self._line0.Initialize()
        self._line0 .InsertNextCell(2)
        self._line0 .InsertCellPoint(0)
        self._line0 .InsertCellPoint(1)
        
        self._line1.Initialize()
        self._line1.InsertNextCell(2)
        self._line1.InsertCellPoint(1)
        self._line1.InsertCellPoint(2)
        self.UpdateDisplay()
    def AddActor(self,ren=None):
        ren.AddActor(self.__actor0)
        ren.AddActor(self.__actor1)
    def UpdateDisplay(self):
        if self.__DisplayType==0:
            self.__mapper0.SetInputConnection(self.__tube0.GetOutputPort())
            self.__mapper1.SetInputConnection(self.__tube1.GetOutputPort())
        else:
            self.__mapper0.SetInputData(self.__pd0)
            self.__mapper1.SetInputData(self.__pd1)
            
        self.__tube0.SetRadius(self.__Radius)
        self.__tube0.SetNumberOfSides(self._Resolution)
        self.__tube1.SetRadius(self.__Radius)
        self.__tube1.SetNumberOfSides(self._Resolution)
        
        self.__actor0.GetProperty().SetOpacity(self.__DisplayOpacity)
        self.__actor0.GetProperty().SetColor(self.__DisplayColor0)
        self.__actor1.GetProperty().SetOpacity(self.__DisplayOpacity)
        self.__actor1.GetProperty().SetColor(self.__DisplayColor1)
    def SetResolution(self,re=20):     
        self._Resolution=re
        self.UpdateDisplay()
    def SetRadius(self,ra=0.2):
        self.__Radius=ra
        self.UpdateDisplay()
    def SetOpacity(self,opac=0.9):
        self.__DisplayOpacity=opac
        self.__actor0.GetProperty().SetOpacity(self.__DisplayOpacity)
        self.__actor1.GetProperty().SetOpacity(self.__DisplayOpacity)
    def SetColor0(self,color=[0,0,0]):
        self.__DisplayColor0=color
        self.__actor0.GetProperty().SetColor(self.__DisplayColor0)
    def GetColor(self):
        return self.__DisplayColor0
    def SetColor1(self,color=[1,1,1]):
        self.__DisplayColor1=color
        self.__actor1.GetProperty().SetColor(self.__DisplayColor1)
    def SetColors(self,color0=[0,0,0],color1=[1,1,1]):
        self.__DisplayColor0=color0
        self.__DisplayColor1=color1
        self.UpdateDisplay()
    def SetDisplayType(self,displaytype=0):
        self.__DisplayType=displaytype
        self.UpdateDisplay()
    def SetVisibility(self,vis):
        self.__actor0.SetVisibility(vis)
        self.__actor1.SetVisibility(vis)
class ViewPeriodBond(ViewBond):
    def __init__(self,parentstructure,atomstart=np.array([[0,0,0],[0,1,0]]),atomend=np.array([[10,10,10],[0.2,1,1]])):
        super(ViewPeriodBond,self).__init__(parentstructure,atomstart, atomend)
        self._PeriodDisplayPos=[]
    def _CaluPeriodDisplayPos(self):
        del self._PeriodDisplayPos[:]
        expandcell=self._parent.GetExpandCell()
        endpos=(expandcell[1]+1.0)
        startpos=expandcell[0]*1.0
        pos=copy(self._StartAtomPos)
        pos=pos+startpos
        while pos[0]>=startpos[0] and pos[0]<=endpos[0]+0.0001:
            pos[1]=self._StartAtomPos[1]+startpos[1]
            while pos[1]>=startpos[1] and pos[1]<=endpos[1]+0.0001:
                pos[2]=self._StartAtomPos[2]+startpos[2]
                while pos[2]>=startpos[2] and pos[2]<=endpos[2]+0.0001:
                    self._PeriodDisplayPos.append(copy(pos))
                    pos[2]+=1.0
                pos[1]+=1.0
            pos[0]+=1.0
    def _CreatePeriodDisplayData(self):
        self._Points.Initialize()     
        self._line0.Initialize()
        self._line1.Initialize()
        for index,startpos in enumerate(self._PeriodDisplayPos):
            endpos=startpos+self._DirectVector
            midposition=startpos+self._DirectVector*0.5
            cartstartpos=self._parent.FracPosToCartPos(startpos)
            cartmidposition=self._parent.FracPosToCartPos(midposition)
            cartendpos=self._parent.FracPosToCartPos(endpos)
            
            self._Points.InsertPoint(index*3+0,cartstartpos)
            self._Points.InsertPoint(index*3+1,cartmidposition)
            self._Points.InsertPoint(index*3+2,cartendpos)

            self._line0 .InsertNextCell(2)
            self._line0 .InsertCellPoint(index*3+0)
            self._line0 .InsertCellPoint(index*3+1)

            self._line1.InsertNextCell(2)
            self._line1.InsertCellPoint(index*3+1)
            self._line1.InsertCellPoint(index*3+2)
    def SetEndPoint(self,startpos,endpos):
        self._StartAtomPos=startpos
        self._DirectVector=endpos-startpos
        self._CaluPeriodDisplayPos()
        self._CreatePeriodDisplayData()
    def UpdateDisplay(self):

        super(ViewPeriodBond,self).UpdateDisplay()
        
class ViewPolyhedra(object):    
    def __init__(self,parentstruc,DataModelPoly=None):
        self._parent=parentstruc
        self._DataPolyhedra=DataModelPoly
        self._CentreSiteId=0
        self.__Label=None
        self.__CentreSite=None
        self.__VertexSites=[]
        self.__VertexCartPos=[]
        self.__ConvexHull=None
        self._Actor=vtk.vtkActor()
        self._Mapper=vtk.vtkDataSetMapper()
        self.__ugrid=vtk.vtkUnstructuredGrid()
        self._Mapper.SetInputData(self.__ugrid)
        self._Actor.SetMapper(self._Mapper)
       # self._Actor.GetProperty().BackfaceCullingOn()
        self.__pts=vtk.vtkPoints()
        self._parent.GetRender().AddActor(self._Actor)
        if DataModelPoly:
            self.SetPolyhedraData(DataModelPoly)
    def SetColor(self,color):
        self._Actor.GetProperty().SetColor(color)
    def GetColor(self):
        return self._Actor.GetProperty().GetColor()
    def SetOpacity(self,opac):
        self._Actor.GetProperty().SetOpacity(opac)
        #self.__Actor.GetProperty().BackfaceCullingOff()
    def GetOpacity(self):
        self._Actor.GetProperty().GetOpacity()
    def SetPolyhedraData(self,DataModelPoly):
        self._DataPolyhedra=DataModelPoly
        self.CreateHull()
    def get_vertex_cart_pos(self):
        return self.__VertexCartPos
    def set_vertex_cart_pos(self, value):
        self.__VertexCartPos = value
    def del_vertex_cart_pos(self):
        del self.__VertexCartPos

    def get_label(self):
        return self.__Label
    def set_label(self, value):
        self.__Label = value
    def del_label(self):
        del self.__Label
        
    def get_vertex_sites(self):
        return self.__VertexSites
    def set_vertex_sites(self, value):
        self.__VertexSites = value
    def del_vertex_sites(self):
        del self.__VertexSites

    def get_centre_site(self):
        return self.__CentreSite
    def set_centre_site(self, value):
        self.__CentreSite = value
    def del_centre_site(self):
        del self.__CentreSite

    def CreateHull(self):

        centrepos=self._DataPolyhedra.CentreSite.GetPosition()
        for dtpos in self._DataPolyhedra.VertexVectors:
            pos=centrepos+dtpos
            self.__VertexCartPos.append(self._parent.FracPosToCartPos(pos))
#        print(self.__VertexCartPos)
        if len(self.__VertexCartPos)<3:
            return
        self.__ConvexHull=ConvexHull(self.__VertexCartPos,qhull_options='QJ')
        for pt  in self.__VertexCartPos:
            self.__pts.InsertNextPoint(pt.tolist())
        self.__ugrid.SetPoints(self.__pts)
        for ctri in self.__ConvexHull.simplices:
            tri=vtk.vtkTriangle()
            tri.GetPointIds().SetId(0,ctri[0])
            tri.GetPointIds().SetId(1,ctri[1])
            tri.GetPointIds().SetId(2,ctri[2])
            self.__ugrid.InsertNextCell(tri.GetCellType(), tri.GetPointIds())
    def AddActor(self,ren=None):
        ren.AddActor(self._Actor)
    def SetVisibility(self,value):
        if value:
            self._Actor.VisibilityOn()
        else:
            self._Actor.VisibilityOff()
    def Hidden(self):
        self._Actor.VisibilityOff()
    def Show(self):
        self._Actor.VisibilityOn()
    Label = property(get_label, set_label, del_label, "Label's docstring")
    CentreSite = property(get_centre_site, set_centre_site, del_centre_site, "CentreSite's docstring")
    VertexSites = property(get_vertex_sites, set_vertex_sites, del_vertex_sites, "VertexSites's docstring")
    VertexCartPos = property(get_vertex_cart_pos, set_vertex_cart_pos, del_vertex_cart_pos, "VertexCartPos's docstring")
    
class ViewPeriodPolyhedra(ViewPolyhedra):
    def __init__(self,parent,DataModelPoly=None):
        self._PeriodDisplayPos=[]
        self._pts={}
        self._ugrids={}
        self._ugriddata=vtk.vtkAppendPolyData()
        super(ViewPeriodPolyhedra,self).__init__(parent,DataModelPoly)
        self._Mapper=vtk.vtkPolyDataMapper()
    def _CaluPeriodDisplayPos(self):
        centrepos=self._DataPolyhedra.CentreSite.GetPosition()
        expandcell=self._parent.GetExpandCell()
        self._PeriodDisplayPos=[]
        endpos=(expandcell[1]+1.0)
        startpos=expandcell[0]*1.0
        pos=copy(centrepos)
        pos=pos+startpos
        while pos[0]>=startpos[0] and pos[0]<=endpos[0]+0.0001:
            pos[1]=centrepos[1]+startpos[1]
            while pos[1]>=startpos[1] and pos[1]<=endpos[1]+0.0001:
                pos[2]=centrepos[2]+startpos[2]
                while pos[2]>=startpos[2] and pos[2]<=endpos[2]+0.0001:
                    self._PeriodDisplayPos.append(copy(pos))
                    pos[2]+=1.0
                pos[1]+=1.0
            pos[0]+=1.0 
    def CreateHull(self):
        self._CaluPeriodDisplayPos()
        VertexCartPos=[]
        self._polydata={}
        for index,centrepos in enumerate(self._PeriodDisplayPos):
            del VertexCartPos[:]
            self._pts[index]=vtk.vtkPoints()
            self._ugrids[index]=vtk.vtkUnstructuredGrid()
            self._polydata[index]=vtk.vtkPolyData()
            for dtpos in self._DataPolyhedra.VertexVectors:
                pos=centrepos+dtpos
                VertexCartPos.append(self._parent.FracPosToCartPos(pos))
    #        print(self.__VertexCartPos)
            if len(VertexCartPos)<4:
                return
            self.__ConvexHull=ConvexHull(VertexCartPos,qhull_options='QJ')
            for pt  in VertexCartPos:
                self._pts[index].InsertNextPoint(pt.tolist())
            self._ugrids[index].SetPoints(self._pts[index])
            for ctri in self.__ConvexHull.simplices:
                tri=vtk.vtkTriangle()
                tri.GetPointIds().SetId(0,ctri[0])
                tri.GetPointIds().SetId(1,ctri[1])
                tri.GetPointIds().SetId(2,ctri[2])
                self._ugrids[index].InsertNextCell(tri.GetCellType(), tri.GetPointIds())
            self._polydata[index].SetPoints(self._pts[index])
            self._polydata[index].SetPolys(self._ugrids[index].GetCells())
            self._ugriddata.AddInputData(self._polydata[index])
        self._Mapper.SetInputConnection(self._ugriddata.GetOutputPort())
        
        
class ViewStructure(object):
    def __init__(self,ren=None,renwin=None,iren=None):
        self._ViewSites=[]
        self._ViewBonds={}
        self._ViewPolyhedras={}
        self._ren=ren
        self._renderWindow=renwin
        self._iren=iren
        self._struct=None
        self._vmw=None;

        self._OutlineActor= vtk.vtkActor()
    def Render(self): 
        self._renderWindow.Render()  
    def GetInactRen(self):
        return self._iren
    def SetInactRen(self,iren):
        self._iren=iren
    def GetRender(self):
        return self._ren
    def SetRender(self,ren=None):  
        self._ren=ren
    def GetRenderWindow(self):
        return self._renderWindow
    def SetRenderWindow(self,renwin=None):
        self._renderWindow=renwin
    def AddViewSite(self,ste=None):
        vst=ViewAtomSite(self,ste)
        self._ViewSites.append(vst)
    def GetViewSites(self):
        return self._ViewSites
    def AddBond(self,bond=None):
        self._ViewBonds.append(bond)
    def CreateDisplay(self):
        self.CreateOutline()
      #  self._iren.SetEnv(self._ren, self._renderWindow)
        self.CreateAxes()
        
        self.CreateBonds()
        self._ren.GetActiveCamera().ParallelProjectionOn()
    def FracPosToCartPos(self,FracPos=None):
        ABC=self._struct.GetABC()
        pos=ABC[0]*FracPos[0]+ABC[1]*FracPos[1]+ABC[2]*FracPos[2]
        return pos
    def CreateBonds(self):
        
        for steId,bonds in self._struct.Bonds.items():
            if  len(bonds)>0:
                self._ViewBonds[steId]=[]
            for bond in bonds:
                bd=ViewPeriodBond(self)
                pos0=bond.Positions[0]
                pos1=bond.Positions[1]
                color0=[i/255.0 for i in self._ViewSites[bond.Site0.GetSiteId()].GetColor()]
                color1=[i/255.0 for i in self._ViewSites[bond.Site1.GetSiteId()].GetColor()]
                bd.SetColors(color0, color1)
                bd.SetEndPoint(pos0,pos1)
                bd.AddActor(self._ren)
                self._ViewBonds[steId].append(bd)
            vpoly=ViewPeriodPolyhedra(self,self._struct.Polyhedras[steId])
            color=[i/255.0 for i in self._ViewSites[steId].GetColor()]
            vpoly.SetColor(color)
            vpoly.CentreSite=self._ViewSites[steId]
            self._ViewPolyhedras[steId]=vpoly
    def GetViewPolyhedras(self):
        return self._ViewPolyhedras
    def GetBonds(self):
        return self._ViewBonds
    def GetStructure(self,stru=None):
        self._struct=stru
        if stru is None:
            print("this is a error structure!!!")
        else:
            ABC=self._struct.GetABC()
            for st in stru.GetUSites():
                self.AddViewSite(st)
            self.CreateDisplay()
    def CreateAxes(self):
        axes=vtk.vtkAxesActor()
        
        vmatrix = vtk.vtkMatrix4x4()
        ABC=self._struct.GetABC().copy()
        lena=math.sqrt(np.sum(ABC[0]**2))
        lenb=math.sqrt(np.sum(ABC[1]**2))
        lenc=math.sqrt(np.sum(ABC[2]**2))
        ABC[0]=ABC[0]/lena
        ABC[1]=ABC[1]/lenb
        ABC[2]=ABC[2]/lenc
        vmatrix.DeepCopy((ABC[0][0],ABC[1][0],ABC[2][0],0,
                                             ABC[0][1],ABC[1][1],ABC[2][1],0,
                                             ABC[0][2],ABC[1][2],ABC[2][2],0,
                                            0,0,0,1))
        transform=vtk.vtkTransform()
        transform.SetMatrix(vmatrix)
        axes.SetUserTransform(transform)
        axes.SetXAxisLabelText('a')
        axes.SetYAxisLabelText('b')
        axes.SetZAxisLabelText('c')

        self._vmw=vtk.vtkOrientationMarkerWidget()
        self._vmw.SetOutlineColor(1,1,1)
        self._vmw.SetInteractor(self._iren)
        self._vmw.SetOrientationMarker(axes)
        self._vmw.SetEnabled(1)
    def ShowAxes(self):
        self._vmw.SetEnabled(1)
    def HideAxes(self):
        self._vmw.SetEnabled(0)
    def HideOutline(self):
        self._OutlineActor.SetVisibility(0)
    def ShowOutline(self):
        self._OutlineActor.SetVisibility(1)
    def CreateOutline(self):
        ABC=self._struct.GetABC()
        pointset=[[0,0,0],ABC[0].tolist(),(ABC[0]+ABC[1]).tolist(),ABC[1].tolist(),[0,0,0]]
        pointset=pointset+[ABC[2].tolist(),(ABC[0]+ABC[2]).tolist(),(ABC[0]+ABC[1]+ABC[2]).tolist(),(ABC[1]+ABC[2]).tolist(),ABC[2].tolist()]
        pointset=pointset+[(ABC[0]+ABC[2]).tolist(),ABC[0].tolist(),(ABC[0]+ABC[1]).tolist(),(ABC[0]+ABC[1]+ABC[2]).tolist()]
        pointset=pointset+[(ABC[1]+ABC[2]).tolist(),ABC[1].tolist()]
        points=vtk.vtkPoints()
        for  point in pointset:
            points.InsertNextPoint(point)
        line2 = vtk.vtkLineSource()
        line2.SetPoints( points ) 
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(line2.GetOutputPort())
        self._OutlineActor.SetMapper(mapper)
        self._OutlineActor.GetProperty().SetOpacity(1.0)
        self._OutlineActor.GetProperty().SetColor([0,0,0])
        self._ren.AddActor(self._OutlineActor)    
    def SetOutlineColor(self, color=[1, 1, 1]):
        self._OutlineActor.GetProperty().SetColor(color)
        
class ViewPeriodStructure(ViewStructure):
    def __init__(self,ren=None,renwin=None,iren=None):
        super(ViewPeriodStructure,self).__init__(ren, renwin, iren)
        self._ExpandCell=np.array([[0,0,0],[0,0,0]])
        
    def AddViewSite(self,ste=None):
        vst=ViewAtomPeriodSite(self,ste)
        self._ViewSites.append(vst)
    def GetExpandCell(self):
        return self._ExpandCell
    def SetExpandCell(self,expandcell):
        self._ExpandCell=expandcell
    


class CustomerRenderWindowInteractor(vtk.vtkGenericRenderWindowInteractor):
    def __init__(self,parent=None):
        self.parent=parent
        self.__Rotating=0
        self.__Panning=0 
        self.__Zooming=0
        self.__renWin=None
        self._ren=None
        self.__ABC=[[1,0,0],[0,1,0],[0,0,1]]

    def SetABC(self,ABC=None): 
        self.__ABC=ABC
    def register(self):
#       self.AddObserver("LeftButtonPressEvent", self.__LeftButtonPressEvent)
#       self.AddObserver("LeftButtonReleaseEvent", self.__LeftButtonReleaseEvent)
#       self.AddObserver("MiddleButtonPressEvent", self.__MiddleButtonPressEvent)
#       self.AddObserver("MiddleButtonReleaseEvent", self.__MiddleButtonReleaseEvent)
#       self.AddObserver("RightButtonPressEvent", self.__RightButtonPressEvent)
#       self.AddObserver("RightButtonReleaseEvent", self.__RightButtonReleaseEvent)
        self.AddObserver("MouseMoveEvent", self.MouseMove)
        self.AddObserver("KeyPressEvent", self.Keypress)
    def SetEnv(self,ren=None,renwin=None):
        self._ren=ren
        self.__renWin=renwin
        self.SetRenderWindow(self.__renWin)
        self._ren.GetActiveCamera().ParallelProjectionOn()
        self.register()
        self.Initialize()
        self._ren.ResetCamera()
    def __LeftButtonPressEvent(self, obj,event):
        if self.__Rotating == 1:
            self.__Rotating=0
        else:
            self.__Rotating=1
    def __LeftButtonReleaseEvent(self, obj,event):
        self.__Rotating=0
    def __MiddleButtonPressEvent(self,obj,event):
            if self.__Panning == 1:
                self.__Panning=0
            else:
                self.__Panning=1
    def __MiddleButtonReleaseEvent(self,obj,event):
            self.__Panning=0
    def __RightButtonPressEvent(self,obj,event):
        if self.__Zooming == 1:
            self.__Zooming=0
        else:
            self.__Zooming=1
    def __RightButtonReleaseEvent(self,obj,event):
        self.__Zooming=1
# General high-level logic
    def MouseMove(self, obj,event):

        lastXYpos = obj.GetLastEventPosition()
        lastX = lastXYpos[0]
        lastY = lastXYpos[1]

        xypos =obj.GetEventPosition()
        x = xypos[0]
        y = xypos[1]

        center = self.__renWin.GetSize()
        centerX = center[0]/2.0
        centerY = center[1]/2.0

        if self.__Rotating:
            self.Rotate(self._ren, self._ren.GetActiveCamera(), x, y, lastX, lastY,
                   centerX, centerY)
        elif self.__Panning:
            self.Pan(self._ren, self._ren.GetActiveCamera(), x, y, lastX, lastY, centerX,
                centerY)
        elif self.__Zooming:
            self.Dolly(self._ren, self._ren.GetActiveCamera(), x, y, lastX, lastY,
                  centerX, centerY)
    def Keypress(self, obj,event):
        key = obj.GetKeySym()

        if key == "w":
            self.Wireframe()
        elif key =="s":
            self.Surface()
        elif key =="r":
            self._ren.ResetCamera()
        elif key =="a":
            self.Viewaaxis()
        elif key =="b":
            self.Viewbaxis()
        elif key =="c":
            self.Viewcaxis()
        elif key =="A":
            self.Viewaaxis(-1.0)
        elif key =="B":
            self.Viewbaxis(-1.0)
        elif key =="C":
            self.Viewcaxis(-1.0)

# Routines that translate the events into camera motions.
    def Viewaaxis(self,axdir=1.0):
            camera=self._ren.GetActiveCamera()
            focal=(self.__ABC[0]+self.__ABC[1]+self.__ABC[2])*0.5
            vp=np.cross(self.__ABC[0],self.__ABC[1])
            camera.SetFocalPoint(focal)
            camera.SetPosition((self.__ABC[0]*5*axdir+focal))
            camera.SetViewUp(vp)
            self._ren.ResetCamera()
            self.__renWin.Render()
    def Viewbaxis(self,axdir=1.0):
            camera=self._ren.GetActiveCamera()
            focal=(self.__ABC[0]+self.__ABC[1]+self.__ABC[2])*0.5
            camera.SetFocalPoint(focal)
            camera.SetPosition((self.__ABC[1]*5*axdir+focal))
            camera.SetViewUp(np.cross(self.__ABC[1],self.__ABC[2]))
            self._ren.ResetCamera()
            self.__renWin.Render()
    def Viewcaxis(self,axdir=1.0):
            camera=self._ren.GetActiveCamera()
            focal=(self.__ABC[0]+self.__ABC[1]+self.__ABC[2])*0.5
            camera.SetFocalPoint(focal)
            camera.SetPosition((self.__ABC[2]*5*axdir+focal))
            camera.SetViewUp(np.cross(self.__ABC[2],self.__ABC[0]))
            self._ren.ResetCamera()
            self.__renWin.Render()
# This one is associated with the left mouse button. It translates x
# and y relative motions into camera azimuth and elevation commands.
#     def Rotate(self,renderer, camera, x, y, lastX, lastY, centerX, centerY):
#         camera.Azimuth(lastX-x)
#         camera.Elevation(lastY-y)
#         camera.OrthogonalizeViewUp()
#         self.__renWin.Render()
# Pan translates x-y motion into translation of the focal point and
# position.
    def Pan(self,renderer, camera, x, y, lastX, lastY, centerX, centerY):
        FPoint = camera.GetFocalPoint()
        FPoint0 = FPoint[0]
        FPoint1 = FPoint[1]
        FPoint2 = FPoint[2]

        PPoint = camera.GetPosition()
        PPoint0 = PPoint[0]
        PPoint1 = PPoint[1]
        PPoint2 = PPoint[2]

        renderer.SetWorldPoint(FPoint0, FPoint1, FPoint2, 1.0)
        renderer.WorldToDisplay()
        DPoint = renderer.GetDisplayPoint()
        focalDepth = DPoint[2]

        APoint0 = centerX+(x-lastX)
        APoint1 = centerY+(y-lastY)

        renderer.SetDisplayPoint(APoint0, APoint1, focalDepth)
        renderer.DisplayToWorld()
        RPoint = renderer.GetWorldPoint()
        RPoint0 = RPoint[0]
        RPoint1 = RPoint[1]
        RPoint2 = RPoint[2]
        RPoint3 = RPoint[3]

        if RPoint3 != 0.0:
            RPoint0 = RPoint0/RPoint3
            RPoint1 = RPoint1/RPoint3
            RPoint2 = RPoint2/RPoint3

        camera.SetFocalPoint( (FPoint0-RPoint0)/2.0 + FPoint0,
                              (FPoint1-RPoint1)/2.0 + FPoint1,
                              (FPoint2-RPoint2)/2.0 + FPoint2)
        camera.SetPosition( (FPoint0-RPoint0)/2.0 + PPoint0,
                            (FPoint1-RPoint1)/2.0 + PPoint1,
                            (FPoint2-RPoint2)/2.0 + PPoint2)
        self.__renWin.Render()


# Dolly converts y-motion into a camera dolly commands.
    def Dolly(self,renderer, camera, x, y, lastX, lastY, centerX, centerY):
        dollyFactor = pow(1.02,(0.5*(y-lastY)))
        if camera.GetParallelProjection():
            parallelScale = camera.GetParallelScale()*dollyFactor
            camera.SetParallelScale(parallelScale)
        else:
            camera.Dolly(dollyFactor)
            renderer.ResetCameraClippingRange()

        self.__renWin.Render()

# Wireframe sets the representation of all actors to wireframe.
    def Wireframe(self):
        actors = self._ren.GetActors()
        actors.InitTraversal()
        actor = actors.GetNextItem()
        while actor:
            actor.GetProperty().SetRepresentationToWireframe()
            actor = actors.GetNextItem()
        self.__renWin.Render()

# Surface sets the representation of all actors to surface.
    def Surface(self):
        actors = self._ren.GetActors()
        actors.InitTraversal()
        actor = actors.GetNextItem()
        while actor:
            actor.GetProperty().SetRepresentationToSurface()
            actor = actors.GetNextItem()
        self.__renWin.Render()
if __name__ == "__main__":
    pass
#     stu=ViewPeriodStructure()
#     vb=ViewPeriodBond(stu)
#     vb.SetEndPoint(np.array([0,0,0]), np.array([0,0.2,0.2]))
