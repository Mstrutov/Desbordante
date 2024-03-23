from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import desbordante

class bcolors:
    HEADER = '\033[95m'
    WARNING = '\033[93m'
    ENDC = '\033[0m'

graph_name = 'channels_graph'
gfd_name = 'entertainment_viewer_gfd'

GRAPH = Path(f'examples/verifying_gfd/datasets/graphs/{graph_name}.dot')
GFD = Path(f'examples/verifying_gfd/datasets/gfds/{gfd_name}.dot')

GRAPH_IMAGE = Path(f'examples/verifying_gfd/figures/graphs/{graph_name}.png')
GFD_IMAGE = Path(f'examples/verifying_gfd/figures/gfds/{gfd_name}.png')

GRAPH_INFO = 'Figure provides an example of a graph dependency and graph. The vertices of the graph have '\
             'labels of C (Channel) or U (User). Depending on the label, the vertex has its own set of attributes. '\
             'In this example, all vertices have a single element list of attributes. At the vertices '\
             'labeled C it consists of the element "topic", and at the vertices labeled U - "age_group". '\
             'On the figure the specific values of these attributes are specified next to the vertices.\n'

GFD_INFO = 'Dependency means that if a user is signed on a channel whose topic is entertainment, he must be a kid.\n'

INFO = 'Let\'s check if this dependency holds.\n'

RESULTS = 'The test found that the constructed dependency is not satisfied because in the '\
          'graph there is an example in which a teenager subscribes to the entertainment channel.\n'

EXIT = f'{bcolors.WARNING}Close the image window to finish.{bcolors.ENDC}'

def execute_algo(algo):
    algo.load_data(gfd=[GFD], graph=GRAPH)
    algo.execute()
    print(f'{bcolors.HEADER}Desbordante > {bcolors.ENDC}', end='')
    if (len(algo.get_gfds()) == 1):
        print('GFD holds.\n')
    else:
        print('GFD does not hold.\n')

def show_example():
    _, axarr = plt.subplots(1,2,figsize=(12,5), gridspec_kw={'width_ratios': [7, 1], 'wspace': 0.5}) 
    axarr[0].set_axis_off()
    axarr[0].set_title(r'$Graph$')
    axarr[0].imshow(mpimg.imread(GRAPH_IMAGE))
    axarr[1].set_axis_off()
    axarr[1].set_title(r'$GFD$')
    axarr[1].imshow(mpimg.imread(GFD_IMAGE))
    axarr[1].text(x=-120, y=250, s=r'$\{0.topic=entertainment \rightarrow$', fontsize=14)
    axarr[1].text(x=-105, y=280, s=r'$1.age\_group=kid\}$', fontsize=14)
    plt.show()

print(GRAPH_INFO)
print(GFD_INFO)
print(INFO)
execute_algo(desbordante.gfd_verification.algorithms.EGfdValid())
print(RESULTS)
print(EXIT)

show_example()
