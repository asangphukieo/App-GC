column(width = 12,
       fluidRow(column(width = 12,
                       mainPanel(width = 12, 
                                 div(class="maincontent",
                                     HTML("<h1>App-GC</h1>"),
                                     HTML("<p>An automated data Preprocessing Pipeline for GC- and GCxGC-MS data</p>"),
                                     HTML("<p>Gas chromatography-mass spectrometry (GC-MS) and comprehensive two-dimensional gas chromatography mass spectrometry (GC×GC-MS) are powerful technologies to carry out untargeted metabolomics profiling and biomarker discovery for biological samples. Data preprocessing is the most crucial step in untargeted metabolomics, since it has direct impact on subsequent statistical analysis and biological interpretation. Therefore, we developed R-based pipeline, App-GC that integrates open-source methods to preprocess GC-MS and GCxGC-MS data. The pipeline consists of two modes, 1D mode for GC-MS data and 2D mode for GCxGC-MS data, each of which contains four main steps 1) peak detection, 2) peak export, 3) peak alignment, and 4) peak annotation.</p>"),
                                     img(src="img/workflow.png", height=500, width=1250, align="center")
                                 )
                                
                       )
       )),#end fluidRow1
       fluidRow(
         column(width=4,
                wellPanel(
                  HTML("<h3>Tutorial</h3>"),
                  HTML("<p>for using App-GC</p>")
                )
         ),
         column(width=4,
                wellPanel(
                  HTML("<h3>Download</h3>"),
                  HTML("<p><a href=https://github.com/asangphukieo/App-GC>App-GC</a> and Example data</p>")
                )
         ),
         column(width=4,
                wellPanel(
                  HTML("<h3>About us</h3>"),
                  HTML("<p><a href=http://metsysbio.com>MSB research group</a></p>")
                )
         ))#end fluidRow2
)#end
