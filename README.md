README OveR~L00P 

This tool is intended to help visualize and explore Distal Interacting Regions (DIR) features that connect in the 3D space of the genome a set of Regions Of Interest (ROI). This readme file is ment to aggregate all the descriptions for the different functions offered by the tool.

######### Initial Loading ##########

- When launched, the first window is an Upload-box tha can take (if everythings fine) any tab-separated text file (tested on TSV & txt) that contains at minimum the following columns:
	- roi_id ==> a concatenation of the coordinates of your ROIs in the followinf format : chrX_startCoord_endCoord
	- distal_chr ==> the seqnames of a DIR for the corresponding ROI
	- distal_start ==> the start coordinate of a DIR
	- distal_end ==> the end coordinate of a DIR
	- Condition ==> The MicroC condition/cell line from which the DIR info come from
	- geneSymbol ==> HGNC name corresponding to the ROI (can be any kind of name)
	- Context ==> If the ROI have been obtained in a specific expermiental context, it should be indivated there, so it can be used as a filtering or comparative parameters later. 
	- overlap_* ==> Pre-calculated overlapping genomic features for the DIR (can be ChIP-seq/ATAC, etc) in a boolean mode (TRUE/FALSE for presence/absence of the feature in this DIR)

- This first input table can be loaded or by clicking on the upload box and choosing the path to the file or by drag and drop.

######### Main Page ##########

On this page several informations are already diplayed at launch. 

First you will find the context(s) that were present in your input file that are, by default, selected but that can be unselected by clicking. If the user want to compare contexts, two modes are proposed, the additive mode (ROI present in context A OR B) and an Intersection mode (ROIs present in A AND B). These options allows to explore the relathionship and behaviour of DIR connected to ROIs between contexts.

The main information is tha alluvial style graph displayed at the center of the page. it gives the possible dynamic relationships between different MicroC conditions for the selected contexts. On top of the graph a Global summary of the number of DIR by selected Context is available. 
These dynamic relathionship can be of 7 types if you have 3 MicroC analyzed:
	- Maintained (DIRs present in the 3 MicroC)
	- Lost Late (DIR present in the two first MicroC then lost)
	- Loss/Regain (DIR present in first and last displayed MicroC but not the middle one)
	- Lost Early (DIRs present only in first MicroC)
	- Gained Early (DIR appearing in second MicroC and maintained with third one)
	- Transient (DIR only present in second MicroC)
	- Gained Late (DIRs only present in last MicroC)

For the displayed DIRs dynamics a feature Analysis table is found below the alluvial graph. In this table user will find the percentage of DIR in the corresponding dynamic condition (column name) that possess the feature found in row name. Column can be sorted by decreasing order of feature presence by clicking on the colname. 

On the left side, below the Contexts selector section, a drop-down section allow user to FIX a feature in the DIRs and create a subselection that when apply update the ROI list to filter only the ROI with at least one DIR in any Condition (selected through the check boxes) that possess this fixed feature. User can fix a feature for its TRUE or FALSE value to give more flexibility to the post analysis. 

An additional option for the alluvial representation is the possibility to switch between the Global View and the Individual ROI View, the global view gives an overview of all DIRs for the selected parameters (Contexts, Features filters), and by clicking on a specific ROI in the scrolling lisst (or by inputing the name of your choice), the alluvial view will automatically switch to the Individual view, showing for the specifically selected ROI its own DIR dynamic and feature composition. 
Still on this Individual ROI View, user can click on the title of the graph (geneSymbol (coordinates)) to open a new window with the genecards info about the selected gene. Along with this option a "Jump to UCSC" button allows to directly go to the (euro mirror of) UCSC to have a more detail view of the region. 

On this main page, 4 additional buttons are usable at the top right corner :

	- Export Excel will create an Excel file with as many tab as MicroC conditions. It will give the subset of selected ROIs by the different applied filters, and the info are separated by conditions.
	- Export SVG will download the Alluvial graph and the table of features in SVG format.
	- Venn diagram button is usable to compare the different contexts (when selection of contexts is comprised between 2 and 4 different contexts) and it will produce the filtered Venn diagram for the selected contexts. 
	- GO Analysis button use the gProfiler API with default parameters to test with Over Representation Analysis the enriched pathway for your selected sets of gene/ROI (4 graphs are displayed GO:BP, GO:MF, GO:CC, ReactomePathway), each table being downloadable by rightclick ==> save image

######### Network Page ##########
On the top left corner, you can click on load DEGs button, which allow to add anlother input table of DiffExpr genes or any additional gene sets different from you ROIs. By doing so you will see a new button appear near the Export Excel one. This "NEtwork View" button allows to drow interaction network between the ROIS selected and the secondary input table. 

Secondary input table format : 

	- geneSymbol column
	- Context column
	- DEG_chr 		|
	- DEG_start		| coordinates of the gene
	- DEG_end		|
	- geneType (coding, lnc, etc)
	- distal_chr	|
	- distal_start	| coordinates of the DIR
	- distal_end	| 

once loaded, you can select the context you want to keep by checking/uncheking the check-boxes before proceeding with further analyses.

Then, by clicking on the NEtwork View button (can take few seconds depending on your machine performance and size of the filtered sets of ROIs) a new window will open with the gene interaction network. Primary loaded Network only goves information about the interacting genes withour more precision. By clicking on Color by Dynamics, you will display the Dynamic interaction in which the genes are connnected. See the paper for more detail about direct and indirect interaction definition.
By mousing over the gene Nodes, you will have information about the context in which the genes are found.

3 more buttons are available on this Network page : 

1/Export SVG button will download the Network figure.
2/Analyse GO will perform a Gene Ontology analysis via gProfiler But, since the network genes are far less than only ROIs, we filter the results by pValue and not FDR/Qvqlue, allowing more output but less powerful statistically speaking. So it can give hint, but to take with cautious consideration. 
3/Close button will do what is expected to do by a button called "close"

