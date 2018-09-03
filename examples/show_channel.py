import vtk
from zeo.netstorage import AtomNetwork
from zeo.channel import Channel

#根据AtomNetwork构造展示AtomNetwork的源数据
#点和线的数据
def AtmnetSource(atomnetwork):
    atoms_site = []
    atoms = atomnetwork.atoms
    for i in atoms:
        atoms_site.append(i[1])
    sphere=vtk.vtkSphereSource()
    sphere.SetCenter(atoms_site)
    sphere.SetRadius(0.04)
    return sphere


#根据Channel构造展示Channel的源数据
#点和线的数据且需要处理周期性
def ChannelSource(channel):
    pass





def Display(atomnetwork,channel):
    
    renWin = vtk.vtkRenderWindow()
    renWin.SetSize(600, 600)
    renWin.SetWindowName('Channel view test')

    renderer = vtk.vtkRenderer()
    
    atoms = atomnetwork.atoms
    #atoms = atomnetwork
    for i in range(len(atoms)):
        atom = vtk.vtkSphereSource()
        atom.SetCenter(atoms[i][1])
        atom.SetRadius(0.04)
        
        atomMapper = vtk.vtkPolyDataMapper()
        atomActor = vtk.vtkActor() 
        
        atomMapper.SetInputConnection(atom.GetOutputPort())
        atomActor.SetMapper(atomMapper)
        atomActor.GetProperty().SetColor(0.9804, 0.9216, 0.8431)
        renderer.AddViewProp(atomActor)
    
    initials = channel.nodes
    #test code
    #initials = [[0, [3.190295850833425, 0.8755507257066707, 2.456838230302852], 1.5088341084285433], [1, [1.1083746175913838, 0.19310377416061847, 5.545615109204309], 1.395878894814294], [2, [0.16072723136374067, 3.9351370561729335, 5.356215528814367], 1.3964483726943062], [3, [1.9748328675382072, 4.033670883989947, 2.1370451794264467], 1.235088970054386], [4, [0.9214256508857486, 1.1748713352351365, 1.5565669858442992], 1.1195247343236383], [5, [2.401335408641232, 3.6907907651916476, 0.10418166251447392], 1.2816998243989537], [6, [0.7388113212217622, 0.9332478908561466, 5.324351976395488], 1.3720916114003], [7, [0.780053889480055, 0.7037759112780848, 5.305620009992378], 1.3422207104331099], [8, [0.8062328350612087, 1.030695449900831, 5.433452109901365], 1.3872567955007895], [9, [1.4342284218232422, 3.562010219494333, 0.3122492683885111], 1.4016814054440634], [10, [1.9859518281340114, 3.6748715011934325, 1.2688627931343748], 1.1195246579752263]]
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(initials))
    for j in range(len(initials)):
        points.InsertPoint(initials[j][0],initials[j][1])   
        
    connections = channel.connections
    #test code
    #connections = [[0,1],[0,2],[0,4],[2,4],[2,5],[3,7],[3,9],[5,6],[5,1],[7,10],[8,9]]
    linkedIdList = vtk.vtkIdList()
    #linkedIdList.InsertNextId(len(connections))
    for k in range(len(connections)):
        #linkedIdList.InsertNextId(2)
        linkedIdList.InsertNextId(connections[k][0])
        linkedIdList.InsertNextId(connections[k][1])
    uGrid = vtk.vtkUnstructuredGrid()
    uGrid.InsertNextCell(vtk.VTK_LINE, linkedIdList)
    uGrid.SetPoints(points)
    
    uGridMapper = vtk.vtkDataSetMapper()
    uGridActor = vtk.vtkActor()
    
    uGridMapper.SetInputData(uGrid)
    uGridActor.SetMapper(uGridMapper)
    uGridActor.GetProperty().SetColor(0.3,0.7,0.9)
    renderer.AddViewProp(uGridActor)
    
    renderer.SetBackground(0.0,0.3,0.1)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(30)
    renderer.GetActiveCamera().Elevation(-30)
    renderer.GetActiveCamera().Zoom(1) 
    renderer.ResetCameraClippingRange()

    renWin.AddRenderer(renderer)

    iRen = vtk.vtkRenderWindowInteractor()
    iRen.SetRenderWindow(renWin)
    iRen.Initialize()
    renWin.Render()
    return iRen

if __name__ == '__main__':
    
    radii = {}
    atmnet = AtomNetwork.read_from_CIF("icsd_16713.cif", radii, False, None)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    channels = Channel.findChannels(vornet,1.1,"icsd_16713.net")
    #print(channels[0].nodes)
    #print(channels[0].connections)

    #test code
    #atmnet = [['Li1+', [0.2029999941587448, 0.550000011920929, 0.3399999141693115]], ['Li1+', [0.7970000058412552, 0.550000011920929, 0.1600000262260437]], ['Li1+', [0.7970000058412552, 0.44999998807907104, 0.6600000262260437]], ['Li1+', [0.2029999941587448, 0.44999998807907104, 0.8399999737739563]], ['Li1+', [0.703000009059906, 0.050000011920928955, 0.3399999141693115]], ['Li1+', [0.296999990940094, 0.050000011920928955, 0.1600000262260437]], ['Li1+', [0.296999990940094, 0.949999988079071, 0.6600000262260437]], ['Li1+', [0.703000009059906, 0.949999988079071, 0.8399999737739563]], ['C4+', [0.0, 0.9429999999701977, 0.75]], ['C4+', [0.0, 0.05700000002980232, 0.25]], ['C4+', [0.5, 0.4429999887943268, 0.75]], ['C4+', [0.5, 0.5569999814033508, 0.25]], ['O2-', [0.0, 0.6870000064373016, 0.75]], ['O2-', [0.0, 0.31299999356269836, 0.25]], ['O2-', [0.5, 0.18700000643730164, 0.75]], ['O2-', [0.5, 0.812999963760376, 0.25]], ['O2-', [0.14499999582767487, 0.06699997186660767, 0.8199999928474426]], ['O2-', [0.8550000041723251, 0.06699997186660767, 0.6800000071525574]], ['O2-', [0.8550000041723251, 0.9330000281333923, 0.18000000715255737]], ['O2-', [0.14499999582767487, 0.9330000281333923, 0.3199999928474426]], ['O2-', [0.6449999809265137, 0.5669999718666077, 0.8199999928474426]], ['O2-', [0.35500001907348633, 0.5669999718666077, 0.6800000071525574]], ['O2-', [0.35500001907348633, 0.4330000877380371, 0.18000000715255737]], ['O2-', [0.6449999809265137, 0.4330000877380371, 0.3199999928474426]]]
    #channels=[1,2]
    
    iRen = renderer = Display(atmnet, channels[0])
    iRen.Start()
