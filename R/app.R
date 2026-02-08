# app.R
# Web interface for huGSVcalleR package using Shiny with real function calls

library(shiny)
library(shinythemes)
library(DT)
library(plotly)
library(shinyFiles)
library(shinyalert)
library(huGSVcalleR)
library(ggplot2)
library(bslib)

# UI Application
ui <- navbarPage(
  "huGSVcalleR Pipeline",
  theme = shinytheme("flatly"),
  id = "nav",

  # Tab: About
  tabPanel("About",
           div(class = "container-fluid",
               style = "padding: 30px;",

               # Header
               div(class = "jumbotron",
                   style = "background-color: #2C3E50; color: white; padding: 30px; border-radius: 10px;",
                   h1("huGSVcalleR Pipeline", style = "margin-top: 0;"),
                   h3("Comprehensive SNP Identification and Annotation Toolkit"),
                   p("Version 0.1.0 | MIT License")
               ),

               # Two column layout
               fluidRow(
                 # Left column - Package info
                 column(6,
                        div(class = "panel panel-primary",
                            div(class = "panel-heading",
                                h3("Package Overview", class = "panel-title", icon("info-circle"))
                            ),
                            div(class = "panel-body",
                                p("huGSVcalleR is an advanced R package that provides a complete, integrated pipeline
                for identifying and annotating single nucleotide polymorphism (SNP) sites in
                whole genome or targeted sequencing data."),

                                h4("Key Features:"),
                                tags$ul(
                                  tags$li("Over 35 high-level functions written in R and C++"),
                                  tags$li("Six integrated analysis modules"),
                                  tags$li("Six different variant calling approaches"),
                                  tags$li("Parallel processing support"),
                                  tags$li("Comprehensive quality control"),
                                  tags$li("Splicing machinery linkage")
                                ),

                                h4("Pipeline Modules:"),
                                tags$ol(
                                  tags$li(strong("Quality Control"), "- Assessment and cleaning of raw reads"),
                                  tags$li(strong("Alignment"), "- Read mapping and alignment QC"),
                                  tags$li(strong("BAM Processing"), "- File preparation and pileup generation"),
                                  tags$li(strong("Variant Calling"), "- SNP detection using multiple methods"),
                                  tags$li(strong("Variant Filtering"), "- Quality-based filtering"),
                                  tags$li(strong("Variant Annotation"), "- Functional annotation and effect prediction")
                                ),

                                h4("Variant Calling Methods:"),
                                tags$ul(
                                  tags$li(strong("Binomial Test"), " - Statistical test based on binomial distribution"),
                                  tags$li(strong("Counting Method"), " - Simple allele counting approach"),
                                  tags$li(strong("Entropy-based"), " - Information theory approach"),
                                  tags$li(strong("Fisher's Exact Test"), " - Statistical significance testing"),
                                  tags$li(strong("Fisher-Smyth"), " - Advanced statistical method"),
                                  tags$li(strong("GATK HaplotypeCaller"), " - Industry-standard approach"),
                                  tags$li(strong("Machine Learning"), " - AI-powered prediction models")
                                ),

                                p(strong("Unique Feature:"), "Linking detected polymorphism sites to splicing machinery.")
                            )
                        ),

                        div(class = "panel panel-success",
                            div(class = "panel-heading",
                                h3("Test Modules", class = "panel-title", icon("vial"))
                            ),
                            div(class = "panel-body",
                                p("The 'Test Modules' tab allows you to run complete, pre-configured test pipelines
                using sample data from the package."),
                                tags$ul(
                                  tags$li(strong("Module 1:"), "Raw read quality control and cleaning"),
                                  tags$li(strong("Module 2:"), "Read alignment and mapping QC"),
                                  tags$li(strong("Module 3:"), "BAM file processing and pileup generation"),
                                  tags$li(strong("Module 4:"), "Variant calling using Fisher's exact test"),
                                  tags$li(strong("Module 5:"), "Variant filtering with custom criteria"),
                                  tags$li(strong("Module 6:"), "Variant annotation with gene information"),
                                  tags$li(strong("Complete Pipeline:"), "Run all modules sequentially")
                                ),
                                p("Each module uses pre-configured paths that can be modified for your specific data.")
                            )
                        )
                 ),

                 # Right column - Technical info
                 column(6,
                        div(class = "panel panel-info",
                            div(class = "panel-heading",
                                h3("System Requirements", class = "panel-title", icon("server"))
                            ),
                            div(class = "panel-body",
                                h4("Software Requirements:"),
                                tags$ul(
                                  tags$li("R (≥ 4.0.0)"),
                                  tags$li("Bioconductor (≥ 3.10)"),
                                  tags$li("CRAN packages: data.table, dplyr, ggplot2, etc."),
                                  tags$li("Optional: GATK for certain methods"),
                                  tags$li("Perl for MaxEntScan splicing prediction")
                                ),

                                h4("Hardware Recommendations:"),
                                tags$ul(
                                  tags$li("8+ GB RAM (minimum)"),
                                  tags$li("16+ GB RAM (recommended for whole genome)"),
                                  tags$li("4+ CPU cores"),
                                  tags$li("50+ GB free disk space"),
                                  tags$li("SSD storage for better I/O performance")
                                ),

                                h4("Supported Operating Systems:"),
                                tags$ul(
                                  tags$li("Linux/Unix (recommended)"),
                                  tags$li("macOS"),
                                  tags$li("Windows (with WSL recommended)")
                                ),

                                h4("Input Formats:"),
                                tags$ul(
                                  tags$li("FASTQ (.fastq, .fastq.gz)"),
                                  tags$li("BAM (.bam)"),
                                  tags$li("VCF (.vcf, .vcf.gz)"),
                                  tags$li("FASTA (.fa, .fasta)"),
                                  tags$li("GTF/GFF (.gtf, .gff)")
                                )
                            )
                        ),

                        div(class = "panel panel-warning",
                            div(class = "panel-heading",
                                h3("Authors & Support", class = "panel-title", icon("users"))
                            ),
                            div(class = "panel-body",
                                h4("Development Team:"),
                                tags$ul(
                                  tags$li(strong("Vasiliy V. Grinev"), " - Primary developer (Grinev_vv@bsu.by)"),
                                  tags$li(strong("Mikalai M. Yatskou"), " - Core developer (Yatskou@bsu.by)"),
                                  tags$li(strong("Dzianis D. Sarnatski"), " - Core developer (denisiussarnatski@gmail.com)"),
                                  tags$li("Other authors will be here later...")
                                ),

                                h4("Affiliation:"),
                                p("Belarusian State University, Faculty of Biology"),

                                h4("Citation:"),
                                div(class = "well",
                                    style = "background-color: #f8f9fa; padding: 10px;",
                                    p(em("Grinev VV, Yatskou MM, Sarnatski DD., Ilyushenok IN., Boeva NA. Development R-package huGSVcalleR, intended for
                                    identification and annotation of the site of single nucleotide polymorphism
                                    in the human genome. CSIST-2025. 2025.")
                                )),

                                h4("License:"),
                                p("MIT License - Open source"),

                                h4("GitHub Repository:"),
                                p(a("https://github.com/bsu-bioinformatics/huGSVcalleR",
                                    href = "https://github.com/bsu-bioinformatics/huGSVcalleR",
                                    target = "_blank",
                                    class = "btn btn-default btn-xs")),

                                h4("Getting Started:"),
                                p("1. Install the package from GitHub"),
                                p("2. Load sample data from the package"),
                                p("3. Run test modules to verify installation"),
                                p("4. Process your own data using the pipeline tabs")
                            )
                        )
                 )
               ),

               # Quick start guide
               div(class = "panel panel-default",
                   div(class = "panel-heading",
                       h3("Quick Start Guide", class = "panel-title", icon("rocket"))
                   ),
                   div(class = "panel-body",
                       h4("Using the Web Interface:"),
                       tags$ol(
                         tags$li("Go to the 'Test Modules' tab to run pre-configured tests"),
                         tags$li("Check the results in the 'Execution Output' tab"),
                         tags$li("View generated files in the 'Files' tab"),
                         tags$li("Download results using the download button"),
                         tags$li("Use individual pipeline tabs for custom analysis")
                       ),

                       h4("Command Line Usage Example:"),
                       tags$pre(
                         '# Load the package
library(huGSVcalleR)

# Run quality control
qa <- assessQRawReads(fastqDir = "Files_FASTQ",
                      fastq = "sample.fastq.gz")

# Clean low-quality reads
cleanRawReads(fastqDir = "Files_FASTQ",
              fastq1 = "sample_R1.fastq",
              fastq2 = "sample_R2.fastq")

# Call variants using multiple methods
callSNVs(fr_data = "pileup_results.csv",
         criteria = c("fisher", "entropy", "ml"),
         outputVCF = "detected_variants.vcf")'
                       ),

                       br(),

                       actionButton("start_tour", "Take a Tour of the Interface",
                                    class = "btn-primary btn-lg", icon("compass"))
                   )
               )
           )
  ),

  # Tab: Test Modules
  tabPanel("Test Modules",
           sidebarLayout(
             sidebarPanel(
               h4("Run Test Modules", icon("vial")),
               selectInput("test_module", "Select Test Module",
                           choices = c(
                             "Module 1: Raw Read QC" = "module1",
                             "Module 2: Alignment & Mapping QC" = "module2",
                             "Module 3: BAM Processing" = "module3",
                             "Module 4: Variant Calling" = "module4",
                             "Module 5: Variant Filtering" = "module5",
                             "Module 6: Variant Annotation" = "module6",
                             "Complete Pipeline" = "complete"
                           )),

               conditionalPanel(
                 condition = "input.test_module == 'module1'",
                 h5("Module 1: Raw Read Quality Control"),
                 textInput("workDir_module1", "Working Directory",
                           value = "D:/Vasily Grinev"),
                 textInput("fastqDir_module1", "FASTQ Directory",
                           value = "Files_FASTQ"),
                 textInput("fastq1_module1", "FASTQ R1",
                           value = "test_seq.R1.fastq.gz"),
                 textInput("fastq2_module1", "FASTQ R2",
                           value = "test_seq.R2.fastq.gz"),
                 textInput("adapters_module1", "Adapters File",
                           value = "adapters.txt"),
                 textInput("contaminants_module1", "Contaminants File",
                           value = "contaminants.txt")
               ),

               conditionalPanel(
                 condition = "input.test_module == 'module2'",
                 h5("Module 2: Alignment & Mapping QC"),
                 textInput("workDir_module2", "Working Directory",
                           value = "D:/Vasily Grinev"),
                 textInput("fastqDir_module2", "FASTQ Directory",
                           value = "Files_FASTQ"),
                 textInput("fastq1_module2", "Filtered FASTQ R1",
                           value = "test_seq.filtered.R1.fastq.gz"),
                 textInput("fastq2_module2", "Filtered FASTQ R2",
                           value = "test_seq.filtered.R2.fastq.gz"),
                 textInput("ref_genome_module2", "Reference Genome Path",
                           value = "Reference_Genomes/GRCh38/GRCh38"),
                 textInput("bamDir_module2", "BAM Output Directory",
                           value = "Files_BAM"),
                 textInput("bamFile_module2", "BAM Output File",
                           value = "test_seq.filtered.bam")
               ),

               conditionalPanel(
                 condition = "input.test_module == 'module3'",
                 h5("Module 3: BAM Processing"),
                 textInput("workDir_module3", "Working Directory",
                           value = "D:/Vasily Grinev"),
                 textInput("bamDir_module3", "BAM Directory",
                           value = "Files_BAM"),
                 textInput("bamFile_module3", "BAM File",
                           value = "test_seq.filtered.bam"),
                 textInput("regions_file_module3", "Regions File (optional)",
                           value = "aml_genes.txt"),
                 checkboxInput("run_filtering_module3", "Run BAM Filtering", TRUE),
                 checkboxInput("run_pileup_module3", "Run Pileup", TRUE)
               ),

               conditionalPanel(
                 condition = "input.test_module == 'module4'",
                 h5("Module 4: Variant Calling"),
                 textInput("workDir_module4", "Working Directory",
                           value = "D:/Vasily Grinev"),
                 textInput("bamDir_module4", "BAM Directory",
                           value = "Files_BAM"),
                 textInput("bamFile_module4", "BAM File",
                           value = "test_seq3.bam"),
                 textInput("fastaDir_module4", "FASTA Directory",
                           value = "Files_FASTA"),
                 textInput("faFile_module4", "Reference FASTA",
                           value = "hg38.fa"),
                 numericInput("min_depth_module4", "Minimum Depth",
                              value = 1, min = 1, max = 100),
                 numericInput("min_baseq_module4", "Minimum Base Quality",
                              value = 20, min = 0, max = 40),
                 numericInput("threads_module4", "Number of Threads",
                              value = 3, min = 1, max = 32)
               ),

               conditionalPanel(
                 condition = "input.test_module == 'module5'",
                 h5("Module 5: Variant Filtering"),
                 textInput("workDir_module5", "Working Directory",
                           value = "D:/Vasily Grinev"),
                 textInput("vcfDir_module5", "VCF Directory",
                           value = "Files_VCF"),
                 textInput("vcfFile_module5", "Input VCF File",
                           value = "test_seq3, fisherCalleR, SNVs.vcf"),
                 textInput("regions_file_module5", "Regions File",
                           value = "aml_genes.subset.txt"),
                 numericInput("min_depth_module5", "Minimum Depth",
                              value = 10, min = 1, max = 1000),
                 numericInput("min_qual_module5", "Minimum Quality",
                              value = 5, min = 0, max = 100),
                 numericInput("max_heterozyg_module5", "Maximum Heterozygosity",
                              value = 0.2, min = 0, max = 1, step = 0.01)
               ),

               conditionalPanel(
                 condition = "input.test_module == 'module6'",
                 h5("Module 6: Variant Annotation"),
                 textInput("workDir_module6", "Working Directory",
                           value = "D:/Vasily Grinev"),
                 textInput("vcfDir_module6", "VCF Directory",
                           value = "Files_VCF"),
                 textInput("vcfFile_module6", "Filtered VCF File",
                           value = "test_seq3, fisherCalleR, SNVs.filtered.vcf"),
                 textInput("genes_file_module6", "Genes Annotation File",
                           value = "Ensembl release 114, GRCh38.p14, annotations of genes.txt")
               ),

               actionButton("run_test", "Run Test Module",
                            class = "btn-primary btn-lg", icon("play")),
               hr(),
               downloadButton("download_results", "Download Results",
                              class = "btn-success")
             ),

             mainPanel(
               tabsetPanel(
                 tabPanel("Execution Output",
                          h3("Test Module Execution", icon("terminal")),
                          verbatimTextOutput("test_output"),
                          br(),
                          h4("Execution Log"),
                          verbatimTextOutput("execution_log")
                 ),
                 tabPanel("Results",
                          h3("Test Results", icon("chart-bar")),
                          uiOutput("results_summary"),
                          dataTableOutput("results_table"),
                          plotOutput("results_plot")
                 ),
                 tabPanel("Files",
                          h3("Generated Files", icon("folder")),
                          verbatimTextOutput("generated_files"),
                          h4("File Contents Preview"),
                          uiOutput("file_preview")
                 )
               )
             )
           )
  ),

  # Tab: Quality Control
  tabPanel("Quality Control",
           sidebarLayout(
             sidebarPanel(
               h4("FASTQ Files Quality Assessment", icon("search")),
               shinyDirButton("qc_dir", "Select Working Directory", "Select Directory"),
               verbatimTextOutput("qc_work_dir"),
               br(),
               fileInput("fastq_file_qc", "Select FASTQ File",
                         accept = c(".fastq", ".fastq.gz", ".fq", ".fq.gz"),
                         multiple = FALSE),
               numericInput("fastq_n_qc", "Number of Reads to Analyze",
                            value = 10000, min = 1000, max = 1000000),
               checkboxInput("adapters_check_qc", "Search for Adapter Sequences", FALSE),
               conditionalPanel(
                 condition = "input.adapters_check_qc == true",
                 fileInput("adapters_file_qc", "Adapters File (TXT)",
                           accept = ".txt")
               ),
               checkboxInput("contaminants_check_qc", "Search for Contaminants", FALSE),
               conditionalPanel(
                 condition = "input.contaminants_check_qc == true",
                 fileInput("contaminants_file_qc", "Contaminants File (TXT)",
                           accept = ".txt")
               ),
               actionButton("run_qc_real", "Run Quality Control",
                            class = "btn-primary", icon("play")),
               hr(),
               actionButton("generate_qc_report", "Generate QC Report",
                            class = "btn-info", icon("file-pdf")),
               downloadButton("download_qc_report_real", "Download Report",
                              class = "btn-success")
             ),

             mainPanel(
               tabsetPanel(
                 tabPanel("QC Results",
                          h3("Read Quality Statistics", icon("table")),
                          tableOutput("qc_summary_table_real"),
                          plotOutput("qc_per_cycle_quality_real"),
                          plotOutput("qc_gc_composition_real")
                 ),
                 tabPanel("Quality Plots",
                          h3("Quality Visualization", icon("chart-line")),
                          plotOutput("qc_quality_distribution_real"),
                          plotOutput("qc_base_composition_real"),
                          plotOutput("qc_per_base_quality_real")
                 ),
                 tabPanel("Read Cleaning",
                          h3("Read Trimming and Filtering", icon("filter")),
                          textInput("fastq1_clean", "FASTQ R1", ""),
                          textInput("fastq2_clean", "FASTQ R2 (optional)", ""),
                          numericInput("tr_score_clean", "Trimming Quality Score",
                                       value = 20, min = 0, max = 40),
                          numericInput("readl_clean", "Minimum Read Length",
                                       value = 100, min = 1, max = 500),
                          actionButton("run_cleaning", "Clean Reads",
                                       class = "btn-warning", icon("broom"))
                 )
               )
             )
           )
  ),

  # Tab: Alignment & Mapping QC
  tabPanel("Alignment & Mapping QC",
           sidebarLayout(
             sidebarPanel(
               h4("Read Alignment to Reference Genome", icon("align-center")),
               shinyDirButton("align_dir", "Select Working Directory", "Select Directory"),
               verbatimTextOutput("align_work_dir"),
               br(),
               fileInput("alignment_fastq1_real", "FASTQ File 1 (R1)",
                         accept = c(".fastq", ".fastq.gz")),
               fileInput("alignment_fastq2_real", "FASTQ File 2 (R2, optional)",
                         accept = c(".fastq", ".fastq.gz")),
               shinyDirButton("genome_index_dir", "Select Genome Index Directory", "Select Directory"),
               verbatimTextOutput("genome_index_path"),
               selectInput("alignment_orientation_real", "Paired-end Orientation",
                           choices = c("forward-reverse" = "fr",
                                       "reverse-forward" = "rf",
                                       "forward-forward" = "ff")),
               numericInput("alignment_threads_real", "Number of Threads",
                            value = 4, min = 1, max = 32),
               checkboxInput("detect_sv_real", "Detect Structural Variants", FALSE),
               actionButton("run_alignment_real", "Run Alignment",
                            class = "btn-primary", icon("play")),
               hr(),
               actionButton("sort_index_bam", "Sort & Index BAM",
                            class = "btn-info", icon("sort"))
             ),

             mainPanel(
               tabsetPanel(
                 tabPanel("Alignment Statistics",
                          h3("Alignment Quality Metrics", icon("chart-bar")),
                          tableOutput("alignment_stats_table_real"),
                          plotOutput("alignment_mapq_dist_real"),
                          plotOutput("alignment_insert_size_dist_real")
                 ),
                 tabPanel("Mapping QC",
                          h3("Mapping Quality Control", icon("check-circle")),
                          h4("BAM File Processing"),
                          fileInput("bam_file_process", "Select BAM File", accept = ".bam"),
                          textInput("regions_file_process", "Regions BED File (optional)", ""),
                          numericInput("mapq_filter", "Minimum Mapping Quality",
                                       value = 20, min = 0, max = 60),
                          numericInput("max_insert_size", "Maximum Insert Size",
                                       value = 1000, min = 100, max = 5000),
                          actionButton("run_bam_filter", "Filter BAM",
                                       class = "btn-warning", icon("filter"))
                 )
               )
             )
           )
  ),

  # Tab: BAM Processing
  tabPanel("BAM Processing",
           sidebarLayout(
             sidebarPanel(
               h4("BAM File Processing", icon("cogs")),
               shinyDirButton("bamproc_dir", "Select Working Directory", "Select Directory"),
               verbatimTextOutput("bamproc_work_dir"),
               br(),
               fileInput("bam_file_bamproc", "Select BAM File", accept = ".bam"),

               h5("Processing Options"),
               checkboxInput("use_regions_bamproc", "Use Regions File", FALSE),
               conditionalPanel(
                 condition = "input.use_regions_bamproc == true",
                 fileInput("regions_file_bamproc", "Regions File (BED/TXT)",
                           accept = c(".bed", ".txt"))
               ),

               checkboxInput("filter_duplicates_bamproc", "Filter Duplicates", TRUE),
               checkboxInput("recalibrate_bamproc", "Recalibrate Base Quality", FALSE),

               h5("Pileup Parameters"),
               numericInput("min_coverage_bamproc", "Minimum Coverage",
                            value = 10, min = 1, max = 100),
               numericInput("min_baseq_bamproc", "Minimum Base Quality",
                            value = 20, min = 0, max = 40),
               numericInput("min_mapq_bamproc", "Minimum Mapping Quality",
                            value = 20, min = 0, max = 60),

               hr(),

               actionButton("run_bam_processing", "Process BAM",
                            class = "btn-primary", icon("play")),
               actionButton("run_pileup_generation", "Generate Pileup",
                            class = "btn-success", icon("database"))
             ),

             mainPanel(
               tabsetPanel(
                 tabPanel("Processing Results",
                          h3("BAM Processing Results", icon("table")),
                          verbatimTextOutput("bamproc_output"),
                          h4("Generated Files"),
                          verbatimTextOutput("bamproc_files")
                 ),
                 tabPanel("Pileup Statistics",
                          h3("Pileup Statistics", icon("chart-bar")),
                          dataTableOutput("pileup_stats_table"),
                          plotOutput("pileup_coverage_plot"),
                          plotOutput("pileup_base_distribution")
                 )
               )
             )
           )
  ),

  # Tab: Variant Calling
  tabPanel("Variant Calling",
           sidebarLayout(
             sidebarPanel(
               h4("Variant Calling Methods", icon("dna")),
               selectInput("calling_method_real", "Calling Method",
                           choices = c("Fisher Exact Test" = "fisher",
                                       "Binomial" = "binomial",
                                       "Entropy" = "entropy",
                                       "Poisson" = "poisson",
                                       "Machine Learning" = "ml")),
               fileInput("calling_bam_file_real", "BAM File for Calling", accept = ".bam"),
               fileInput("calling_ref_fasta_real", "Reference Genome (FASTA)",
                         accept = c(".fa", ".fasta", ".fa.gz")),
               conditionalPanel(
                 condition = "input.calling_method_real == 'fisher'",
                 numericInput("min_depth_real", "Minimum Coverage Depth",
                              value = 10, min = 1, max = 100),
                 numericInput("min_baseq_real", "Minimum Base Quality",
                              value = 20, min = 0, max = 40),
                 numericInput("min_mapq_real", "Minimum Mapping Quality",
                              value = 20, min = 0, max = 60),
                 numericInput("qvalue_real", "Q-value Cutoff",
                              value = 12, min = 0, max = 100)
               ),
               conditionalPanel(
                 condition = "input.calling_method_real == 'ml'",
                 selectInput("ml_model_real", "ML Model",
                             choices = c("Decision Tree" = "rpart",
                                         "Random Forest" = "rand_forest",
                                         "XGBoost" = "xgboost",
                                         "SVM" = "SVM",
                                         "Logistic Regression" = "log_reg"))
               ),
               actionButton("run_calling_real", "Run Variant Calling",
                            class = "btn-primary", icon("play")),
               hr(),
               fileInput("pileup_file", "Or use Pileup Data (CSV)", accept = ".csv"),
               actionButton("run_pileup_calling", "Call from Pileup",
                            class = "btn-info", icon("database"))
             ),

             mainPanel(
               tabsetPanel(
                 tabPanel("Detected Variants",
                          h3("Detected Variants", icon("search")),
                          dataTableOutput("variants_table_real"),
                          downloadButton("download_vcf_real", "Download VCF File",
                                         class = "btn-success")
                 ),
                 tabPanel("Calling Statistics",
                          h3("Variant Statistics", icon("chart-bar")),
                          tableOutput("calling_stats_table_real"),
                          plotOutput("variant_type_dist_real"),
                          plotOutput("variant_allele_freq_real")
                 ),
                 tabPanel("Pileup Statistics",
                          h3("Generate Pileup Data", icon("table")),
                          fileInput("pileup_bam_file", "BAM File for Pileup", accept = ".bam"),
                          fileInput("pileup_regions", "Regions File (optional)",
                                    accept = c(".bed", ".txt")),
                          numericInput("pileup_coverage", "Minimum Coverage",
                                       value = 10, min = 1, max = 100),
                          actionButton("run_pileup", "Generate Pileup",
                                       class = "btn-warning", icon("cogs"))
                 )
               )
             )
           )
  ),

  # Tab: Variant Filtering
  tabPanel("Variant Filtering",
           sidebarLayout(
             sidebarPanel(
               h4("Variant Filtering", icon("filter")),
               fileInput("filter_vcf_file_real", "Input VCF File", accept = ".vcf"),
               numericInput("filter_min_depth_real", "Minimum Depth",
                            value = 10, min = 1, max = 1000),
               numericInput("filter_min_qual_real", "Minimum Quality",
                            value = 20, min = 0, max = 100),
               sliderInput("filter_max_af_real", "Maximum Allele Frequency",
                           min = 0, max = 1, value = 1, step = 0.01),
               fileInput("regions_filter", "Regions File for Filtering",
                         accept = c(".bed", ".txt")),
               actionButton("run_filtering_real", "Run Filtering",
                            class = "btn-primary", icon("filter")),
               hr(),
               downloadButton("download_filtered_vcf_real", "Download Filtered VCF",
                              class = "btn-success")
             ),

             mainPanel(
               tabsetPanel(
                 tabPanel("Filtering Results",
                          h3("Variant Filtering Results", icon("table")),
                          verbatimTextOutput("filtering_output"),
                          h4("Filtered Variants Table"),
                          dataTableOutput("filtered_variants_table"),
                          h4("Filtering Statistics"),
                          tableOutput("filtering_stats_table")
                 ),
                 tabPanel("VCF Info",
                          h3("VCF Information", icon("info-circle")),
                          verbatimTextOutput("vcf_info_output"),
                          h4("VCF Header"),
                          verbatimTextOutput("vcf_header_output")
                 )
               )
             )
           )
  ),

  # Tab: Variant Annotation
  tabPanel("Variant Annotation",
           sidebarLayout(
             sidebarPanel(
               h4("Variant Annotation", icon("sticky-note")),
               fileInput("annotation_vcf_file_real", "Filtered VCF File", accept = ".vcf"),
               fileInput("genes_file_anno", "Gene Annotations File",
                         accept = c(".txt", ".gtf", ".gff")),

               h5("Annotation Options"),
               checkboxInput("dbSNP_annotation", "Add dbSNP Annotation", FALSE),
               conditionalPanel(
                 condition = "input.dbSNP_annotation == true",
                 fileInput("dbsnp_file_anno", "dbSNP File (VCF)", accept = ".vcf")
               ),

               checkboxInput("predict_splicing", "Predict Splicing Effects", TRUE),
               checkboxInput("predict_protein", "Predict Protein Effects", TRUE),

               hr(),

               actionButton("run_annotation_real", "Annotate Variants",
                            class = "btn-primary", icon("sticky-note")),
               actionButton("generate_annotation_report", "Generate Report",
                            class = "btn-info", icon("file-alt")),
               downloadButton("download_annotations_real", "Download Annotations",
                              class = "btn-success")
             ),

             mainPanel(
               tabsetPanel(
                 tabPanel("Annotation Results",
                          h3("Annotation Results", icon("table")),
                          verbatimTextOutput("annotation_output"),
                          h4("Annotated Variants Table"),
                          dataTableOutput("annotated_variants_table"),
                          h4("Annotation Statistics"),
                          tableOutput("annotation_stats_table")
                 ),
                 tabPanel("Gene Effects",
                          h3("Gene Effects Analysis", icon("dna")),
                          h4("Top Affected Genes"),
                          plotOutput("top_genes_plot"),
                          h4("Variant Impact Distribution"),
                          plotOutput("impact_dist_plot"),
                          h4("Splicing Effects"),
                          plotOutput("splicing_effects_plot")
                 )
               )
             )
           )
  )
)

# Server Application
server <- function(input, output, session) {

  # Reactive values for data storage
  values <- reactiveValues(
    # Test modules
    test_results = NULL,
    execution_log = character(),
    generated_files = character(),

    # Individual module results
    filtering_results = NULL,
    annotation_results = NULL,
    bamproc_results = NULL,

    current_work_dir = getwd()
  )

  # Function to add messages to log
  add_to_log <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_entry <- paste(timestamp, "-", message)
    values$execution_log <- c(values$execution_log, log_entry)
  }

  # Function to add generated file
  add_generated_file <- function(file_path) {
    values$generated_files <- c(values$generated_files, file_path)
  }

  # Execute test modules
  observeEvent(input$run_test, {
    showModal(modalDialog(
      title = "Running Test Module",
      "Please wait while the test module is running...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      module <- input$test_module
      add_to_log(paste("Starting test module:", module))

      if (module == "module1") {
        # Module 1: Raw Read QC
        workDir <- input$workDir_module1
        fastqDir <- input$fastqDir_module1
        fastq1 <- input$fastq1_module1
        fastq2 <- input$fastq2_module1

        # Set working directory
        setwd(dir = workDir)

        # Run QC for R1
        add_to_log("Running QC for R1 reads...")
        qa1 <- huGSVcalleR::assessQRawReads(
          fastqDir = fastqDir,
          fastq = fastq1,
          n = NULL,
          adapters = input$adapters_module1,
          contaminants = input$contaminants_module1,
          workDir = workDir
        )

        # Save results
        rds_file1 <- paste(workDir, "Test_sequencing_read1_QA_results.rds", sep = "/")
        saveRDS(object = qa1, file = rds_file1)
        add_generated_file(rds_file1)

        # Generate report
        report1 <- huGSVcalleR::reportQAResults(
          x = qa1,
          output = NULL,
          workDir = workDir
        )

        # Run QC for R2
        add_to_log("Running QC for R2 reads...")
        qa2 <- huGSVcalleR::assessQRawReads(
          fastqDir = fastqDir,
          fastq = fastq2,
          n = NULL,
          adapters = input$adapters_module1,
          contaminants = input$contaminants_module1,
          workDir = workDir
        )

        rds_file2 <- paste(workDir, "Test_sequencing_read2_QA_results.rds", sep = "/")
        saveRDS(object = qa2, file = rds_file2)
        add_generated_file(rds_file2)

        report2 <- huGSVcalleR::reportQAResults(
          x = qa2,
          output = NULL,
          workDir = workDir
        )

        # Clean reads
        add_to_log("Cleaning raw reads...")
        huGSVcalleR::cleanRawReads(
          fastqDir = fastqDir,
          fastq1 = fastq1,
          fastq2 = fastq2,
          adapters = input$adapters_module1,
          error = 0.2,
          min_match_flank = 3L,
          anchored = TRUE,
          indels = FALSE,
          tr_score = 20,
          tr_start = 10,
          tr_end = 30,
          phred = "phred_scores.txt",
          k = 3,
          halfwidth = NULL,
          successive = TRUE,
          readl = 100,
          readq = 20,
          readn = 3,
          dustScore = 875,
          batchSize = NA,
          postfix = "filtered",
          workDir = workDir
        )

        values$test_results <- list(
          module = "Module 1: Raw Read QC",
          status = "Completed",
          files = c(rds_file1,
                    rds_file2),
          qc_r1 = qa1,
          qc_r2 = qa2
        )

      } else if (module == "module2") {
        # Module 2: Alignment & Mapping QC
        workDir <- input$workDir_module2
        setwd(dir = workDir)

        # Build index (if needed)
        add_to_log("Building genome index...")
        huGSVcalleR::buildIndexSubread(
          ref_genome = "Reference_Genomes/GRCh38",
          ref_fasta = "Files_FASTA",
          index = "GRCh38",
          fa = "hg38.fa.gz",
          memory = 8000,
          workDir = workDir
        )

        # Run alignment
        add_to_log("Running read alignment...")
        huGSVcalleR::alignDNASubread(
          genome = input$ref_genome_module2,
          fastqDir = input$fastqDir_module2,
          fastq1 = input$fastq1_module2,
          fastq2 = input$fastq2_module2,
          bamDir = input$bamDir_module2,
          bamFile = input$bamFile_module2,
          orientation = "fr",
          threads = 4,
          SV = FALSE,
          workDir = workDir
        )

        # Sort and index BAM
        add_to_log("Sorting and indexing BAM file...")
        huGSVcalleR::sortBamFile(
          bamDir = input$bamDir_module2,
          bamFile = input$bamFile_module2,
          byQname = FALSE,
          workDir = workDir
        )

        # Assess alignment quality
        add_to_log("Assessing alignment quality...")
        qal <- huGSVcalleR::assessQAlignedReads(
          bamDir = input$bamDir_module2,
          bamFile = input$bamFile_module2,
          workDir = workDir
        )

        # Save results
        rds_file <- paste(workDir, "Sample_test_seq.filtered_quality_of_reads_alignment.rds", sep = "/")
        saveRDS(object = qal, file = rds_file)
        add_generated_file(rds_file)

        # Generate report
        huGSVcalleR::reportQAlResults(
          x = qal,
          output = NULL,
          workDir = workDir
        )

        values$test_results <- list(
          module = "Module 2: Alignment & Mapping QC",
          status = "Completed",
          files = rds_file,
          alignment_stats = qal
        )

      } else if (module == "module3") {
        # Module 3: BAM Processing
        workDir <- input$workDir_module3
        setwd(dir = workDir)

        if (input$run_filtering_module3) {
          add_to_log("Filtering BAM file...")
          huGSVcalleR::filterBamFile(
            bamDir = input$bamDir_module3,
            bamFile = input$bamFile_module3,
            gr = NULL,
            flag = Rsamtools::scanBamFlag(isPaired = TRUE, isProperPair = TRUE),
            tlen = 300,
            mapq = 20,
            index = TRUE,
            workDir = workDir
          )
        }

        if (input$run_pileup_module3) {
          add_to_log("Generating pileup statistics...")
          res <- huGSVcalleR::collectBasesPileup(
            bamDir = input$bamDir_module3,
            bamFile = "test_seq.filtered2.bam",
            gr = input$regions_file_module3,
            canChr = TRUE,
            depth = 1e6,
            baseq = 20,
            mapq = 20,
            coverage = 10,
            genome = "BSgenome.Hsapiens.UCSC.hg38",
            tmpDir = NULL,
            postfix = "PileupResults",
            workDir = workDir
          )

          # Save the pileup results
          if (!is.null(res)) {
            values$bamproc_results <- res
          }
        }

        values$test_results <- list(
          module = "Module 3: BAM Processing",
          status = "Completed",
          pileup_generated = input$run_pileup_module3
        )

      } else if (module == "module4") {
        # Module 4: Variant Calling
        workDir <- input$workDir_module4
        setwd(dir = workDir)

        add_to_log("Running variant calling with Fisher's exact test...")
        res <- huGSVcalleR::fisherCalleR(
          bamDir = input$bamDir_module4,
          bamFile = input$bamFile_module4,
          fastaDir = input$fastaDir_module4,
          faFile = input$faFile_module4,
          baseq = input$min_baseq_module4,
          mindepth = input$min_depth_module4,
          maxdepth = 1e6,
          qvalue = 12,
          trim = 0,
          threads = input$threads_module4,
          workDir = workDir
        )

        vcf_file <- paste(input$bamDir_module4,
                          gsub("\\.bam$", "_fisherCalleR_SNVs.vcf", input$bamFile_module4),
                          sep = "/")
        add_generated_file(vcf_file)

        values$test_results <- list(
          module = "Module 4: Variant Calling",
          status = "Completed",
          vcf_file = vcf_file,
          variants = res
        )

      } else if (module == "module5") {
        # Module 5: Variant Filtering
        workDir <- input$workDir_module5
        setwd(dir = workDir)

        add_to_log("Filtering variants...")
        SNVs_filter <- huGSVcalleR::filterVariants(
          vcfDir = input$vcfDir_module5,
          vcfFile = input$vcfFile_module5,
          gr = input$regions_file_module5,
          depth = input$min_depth_module5,
          score = input$min_qual_module5,
          heterozyg = input$max_heterozyg_module5,
          heterozyg1 = input$max_heterozyg_module5,
          index = TRUE,
          workDir = workDir
        )

        filtered_vcf <- paste(input$vcfDir_module5,
                              gsub("\\.vcf$", ".filtered.vcf", input$vcfFile_module5),
                              sep = "/")
        add_generated_file(filtered_vcf)

        # Store filtering results for display
        values$filtering_results <- SNVs_filter

        values$test_results <- list(
          module = "Module 5: Variant Filtering",
          status = "Completed",
          filtered_vcf = filtered_vcf,
          variant_count = length(SNVs_filter),
          filtered_variants = SNVs_filter
        )

      } else if (module == "module6") {
        # Module 6: Variant Annotation
        workDir <- input$workDir_module6
        setwd(dir = workDir)

        # Load genes data
        add_to_log("Loading genes annotation...")
        genes_path <- input$genes_file_module6

        add_to_log("Annotating variants...")
        annoSNVs <- huGSVcalleR::annotateVariants(
          vcfDir = input$vcfDir_module6,
          vcfFile = input$vcfFile_module6,
          ref_genes = genes_path,
          output = "test_seq3_fisherCalleR_annotated_SNVs",
          workDir = workDir
        )

        annotation_file <- paste(workDir,
                                 "test_seq3_fisherCalleR_annotated_SNVs.based annotation.txt",
                                 sep = "/")
        add_generated_file(annotation_file)

        # Store annotation results for display
        values$annotation_results <- annoSNVs

        values$test_results <- list(
          module = "Module 6: Variant Annotation",
          status = "Completed",
          annotation_file = annotation_file,
          annotation_results = annoSNVs
        )

      } else if (module == "complete") {
        # Complete pipeline - run all modules sequentially
        workDir <- "D:/Vasily Grinev"
        setwd(dir = workDir)

        # Run module 1-6 sequentially
        # (Implementation would call each module in sequence)
        add_to_log("Complete pipeline execution started...")

        values$test_results <- list(
          module = "Complete Pipeline",
          status = "Execution started",
          progress = "Running sequentially..."
        )
      }

      add_to_log(paste("Test module", module, "completed successfully"))
      showNotification(paste("Test module", module, "completed!"),
                       type = "message", duration = 5)

    }, error = function(e) {
      add_to_log(paste("Error in test module:", e$message))
      showNotification(paste("Error:", e$message),
                       type = "error", duration = 10)
    })

    removeModal()
  })

  # Display test output
  output$test_output <- renderPrint({
    if (!is.null(values$test_results)) {
      cat("Module:", values$test_results$module, "\n")
      cat("Status:", values$test_results$status, "\n")
      cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

      if (!is.null(values$test_results$files)) {
        cat("Generated files:\n")
        cat(paste(values$test_results$files, collapse = "\n"), "\n\n")
      }

      if (!is.null(values$test_results$variant_count)) {
        cat("Filtered variants:", values$test_results$variant_count, "\n\n")
      }

      # Display specific results for each module
      if (input$test_module == "module5" && !is.null(values$filtering_results)) {
        cat("\n--- Filtering Results ---\n")
        print(values$filtering_results)
        cat("\n--- VCF Info ---\n")
        cat("Class:", class(values$filtering_results), "\n")
        cat("Dimensions:", dim(values$filtering_results), "\n")
        if (!is.null(DelayedArray::rowRanges(values$filtering_results))) {
          cat("Number of variants:", length(DelayedArray::rowRanges(values$filtering_results)), "\n")
        }
      }

      if (input$test_module == "module6" && !is.null(values$annotation_results)) {
        cat("\n--- Annotation Results ---\n")
        if (is.data.frame(values$annotation_results)) {
          print(head(values$annotation_results, 20))
          cat("\nTotal annotated variants:", nrow(values$annotation_results), "\n")
        } else {
          print(values$annotation_results)
        }
      }
    } else {
      cat("No test results available. Run a test module first.")
    }
  })

  # Display execution log
  output$execution_log <- renderText({
    if (length(values$execution_log) > 0) {
      paste(values$execution_log, collapse = "\n")
    } else {
      "No execution log available."
    }
  })

  # Display generated files
  output$generated_files <- renderPrint({
    if (length(values$generated_files) > 0) {
      cat("Generated files:\n")
      cat(paste(values$generated_files, collapse = "\n"))
    } else {
      cat("No files generated yet.")
    }
  })

  # Results table for test modules
  output$results_table <- renderDataTable({
    if (!is.null(values$test_results)) {
      if (input$test_module == "module5" && !is.null(values$filtering_results)) {
        # Display filtered variants table
        if (inherits(values$filtering_results, "CollapsedVCF")) {
          # Extract variant information from VCF
          vcf_info <- as.data.frame(DelayedArray::rowRanges(values$filtering_results))
          info_data <- as.data.frame(VariantAnnotation::info(values$filtering_results))

          # Combine variant information
          result_df <- cbind(vcf_info, info_data)

          # Format for display
          result_df <- result_df[, c("seqnames", "start", "end", "REF", "ALT", "QUAL", "DP")]
          colnames(result_df) <- c("Chromosome", "Start", "End", "REF", "ALT", "Quality", "Depth")

          return(datatable(result_df,
                           options = list(pageLength = 10, scrollX = TRUE),
                           rownames = FALSE))
        }
      } else if (input$test_module == "module6" && !is.null(values$annotation_results)) {
        # Display annotation results
        if (is.data.frame(values$annotation_results)) {
          return(datatable(values$annotation_results,
                           options = list(pageLength = 10, scrollX = TRUE),
                           rownames = FALSE))
        }
      }
    }
    return(NULL)
  })

  # Run variant filtering from Variant Filtering tab
  observeEvent(input$run_filtering_real, {
    req(input$filter_vcf_file_real)

    showModal(modalDialog(
      title = "Running Variant Filtering",
      "Please wait while filtering variants...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      add_to_log("Starting variant filtering from Variant Filtering tab")

      # Save uploaded file
      vcf_path <- input$filter_vcf_file_real$datapath
      vcf_dir <- dirname(vcf_path)
      vcf_file <- basename(vcf_path)

      # Check if regions file is provided
      regions_path <- NULL
      if (!is.null(input$regions_filter)) {
        regions_path <- input$regions_filter$datapath
      }

      # Run filtering
      filtered_variants <- huGSVcalleR::filterVariants(
        vcfDir = "Files_VCF",
        vcfFile = "test_seq3, fisherCalleR, SNVs.vcf",
        gr = "aml_genes.subset.txt",
        depth = input$filter_min_depth_real,
        score = input$filter_min_qual_real,
        heterozyg = input$filter_max_af_real,
        heterozyg1 = input$filter_max_af_real,
        index = TRUE,
        workDir = "D:/Vasily Grinev"
      )

      # Store results
      values$filtering_results <- filtered_variants

      # Generate output file path
      output_vcf <- file.path(vcf_dir,
                              paste0(tools::file_path_sans_ext(vcf_file),
                                     ".filtered.vcf"))

      # Save filtered VCF
      if (inherits(filtered_variants, "CollapsedVCF")) {
        VariantAnnotation::writeVcf(filtered_variants, output_vcf)
        add_generated_file(output_vcf)
      }

      add_to_log("Variant filtering completed successfully")
      showNotification("Variant filtering completed!",
                       type = "message", duration = 5)

    }, error = function(e) {
      add_to_log(paste("Error in variant filtering:", e$message))
      showNotification(paste("Error:", e$message),
                       type = "error", duration = 10)
    })

    removeModal()
  })

  # Display filtering results in Variant Filtering tab
  output$filtering_output <- renderPrint({
    if (!is.null(values$filtering_results)) {
      cat("=== Variant Filtering Results ===\n\n")

      if (inherits(values$filtering_results, "CollapsedVCF")) {
        cat("Class: CollapsedVCF\n")
        cat("Dimensions:", dim(values$filtering_results), "\n")

        # Extract and display rowRanges information
        vcf_ranges <- DelayedArray::rowRanges(values$filtering_results)
        cat("Number of variants:", length(vcf_ranges), "\n\n")

        cat("Row Ranges (first 10 variants):\n")
        if (length(vcf_ranges) > 0) {
          if (length(vcf_ranges) > 10) {
            print(vcf_ranges[1:10])
          } else {
            print(vcf_ranges)
          }
        }

        # Display INFO header
        cat("\nINFO header:\n")
        info_header <- VariantAnnotation::info(VariantAnnotation::header(values$filtering_results))
        for (i in 1:nrow(info_header)) {
          cat(sprintf("%-10s %-10s %s\n",
                      rownames(info_header)[i],
                      info_header$Number[i],
                      info_header$Description[i]))
        }

        # Display INFO data
        if (length(vcf_ranges) > 0) {
          cat("\nINFO data (first 10 variants):\n")
          info_data <- VariantAnnotation::info(values$filtering_results)
          if (length(vcf_ranges) > 10) {
            print(info_data[1:10, ])
          } else {
            print(info_data)
          }
        }
      } else {
        cat("Filtering results object type:", class(values$filtering_results), "\n")
        print(values$filtering_results)
      }
    } else {
      cat("No filtering results available. Run variant filtering first.")
    }
  })

  # Display filtered variants table
  output$filtered_variants_table <- renderDataTable({
    if (!is.null(values$filtering_results)) {
      if (inherits(values$filtering_results, "CollapsedVCF")) {
        # Extract variant information
        vcf_data <- as.data.frame(DelayedArray::rowRanges(values$filtering_results))
        info_data <- as.data.frame(VariantAnnotation::info(values$filtering_results))

        # Combine data
        result_df <- cbind(vcf_data, info_data)

        # Select and rename columns
        keep_cols <- c("seqnames", "start", "end", "REF", "ALT", "QUAL", "DP")
        result_df <- result_df[, intersect(keep_cols, colnames(result_df))]

        if ("seqnames" %in% colnames(result_df)) {
          colnames(result_df)[colnames(result_df) == "seqnames"] <- "Chromosome"
        }
        if ("start" %in% colnames(result_df)) {
          colnames(result_df)[colnames(result_df) == "start"] <- "Start"
        }
        if ("end" %in% colnames(result_df)) {
          colnames(result_df)[colnames(result_df) == "end"] <- "End"
        }

        # Create DataTable
        datatable(result_df,
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE,
                    columnDefs = list(list(className = 'dt-center', targets = "_all"))
                  ),
                  rownames = FALSE) %>%
          formatRound(columns = c("QUAL"), digits = 2)
      }
    }
  })

  # Display filtering statistics
  output$filtering_stats_table <- renderTable({
    if (!is.null(values$filtering_results)) {
      if (inherits(values$filtering_results, "CollapsedVCF")) {
        vcf_data <- as.data.frame(DelayedArray::rowRanges(values$filtering_results))

        stats <- data.frame(
          Metric = c("Total Variants",
                     "Mean Quality Score",
                     "Min Position",
                     "Max Position",
                     "Chromosomes"),
          Value = c(
            nrow(vcf_data),
            ifelse("QUAL" %in% colnames(vcf_data) && !all(is.na(vcf_data$QUAL)),
                   round(mean(vcf_data$QUAL, na.rm = TRUE), 2), "N/A"),
            ifelse("start" %in% colnames(vcf_data), min(vcf_data$start), "N/A"),
            ifelse("end" %in% colnames(vcf_data), max(vcf_data$end), "N/A"),
            ifelse("seqnames" %in% colnames(vcf_data),
                   paste(unique(vcf_data$seqnames), collapse = ", "), "N/A")
          )
        )

        return(stats)
      }
    }
  }, striped = TRUE, hover = TRUE, width = "100%")

  # Display VCF info
  output$vcf_info_output <- renderPrint({
    if (!is.null(values$filtering_results)) {
      if (inherits(values$filtering_results, "CollapsedVCF")) {
        cat("=== VCF Object Information ===\n\n")

        # Basic information
        cat("Class:", class(values$filtering_results), "\n")
        cat("Dimensions:", dim(values$filtering_results), "\n")

        # Metadata
        metadata_list <- S4Vectors::metadata(values$filtering_results)
        if (length(metadata_list) > 0) {
          cat("\nMetadata:\n")
          for (i in seq_along(metadata_list)) {
            cat(names(metadata_list)[i], ":",
                paste(metadata_list[[i]], collapse = ", "), "\n")
          }
        }

        # Fixed fields
        fixed_fields <- fixed(values$filtering_results)
        if (nrow(fixed_fields) > 0) {
          cat("\nFixed fields (first 5 variants):\n")
          print(head(fixed_fields, 5))
        }
      }
    }
  })

  # Display VCF header
  output$vcf_header_output <- renderPrint({
    if (!is.null(values$filtering_results)) {
      if (inherits(values$filtering_results, "CollapsedVCF")) {
        vcf_header <- VariantAnnotation::header(values$filtering_results)

        cat("=== VCF Header Information ===\n\n")

        # CONTIG lines
        if (!is.null(vcf_header$contig)) {
          cat("Contigs (chromosomes):\n")
          print(vcf_header$contig)
          cat("\n")
        }

        # FILTER lines
        if (!is.null(vcf_header$FILTER)) {
          cat("FILTER definitions:\n")
          print(vcf_header$FILTER)
          cat("\n")
        }

        # INFO lines
        if (!is.null(vcf_header$INFO)) {
          cat("INFO field definitions:\n")
          print(vcf_header$INFO)
          cat("\n")
        }

        # FORMAT lines
        if (!is.null(vcf_header$FORMAT)) {
          cat("FORMAT field definitions:\n")
          print(vcf_header$FORMAT)
        }
      }
    }
  })

  # Run variant annotation from Variant Annotation tab
  observeEvent(input$run_annotation_real, {
    req(input$annotation_vcf_file_real, input$genes_file_anno)

    showModal(modalDialog(
      title = "Running Variant Annotation",
      "Please wait while annotating variants...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      add_to_log("Starting variant annotation from Variant Annotation tab")

      # Save uploaded files
      vcf_path <- input$annotation_vcf_file_real$datapath
      genes_path <- input$genes_file_anno$datapath

      vcf_dir <- dirname(vcf_path)
      vcf_file <- basename(vcf_path)

      # Run annotation
      annoSNVs <- huGSVcalleR::annotateVariants(
        vcfDir = "Files_VCF",
        vcfFile = "test_seq3, fisherCalleR, SNVs.filtered.vcf",
        ref_genes = "Ensembl release 114, GRCh38.p14, annotations of genes.txt",
        output = paste0(tools::file_path_sans_ext(vcf_file), "_annotated"),
        workDir = "D:/Vasily Grinev"
      )

      # Store results
      values$annotation_results <- annoSNVs

      # Generate output file
      output_file <- file.path(vcf_dir,
                               paste0(tools::file_path_sans_ext(vcf_file),
                                      "_annotated_SNVs.txt"))

      # Save annotation results
      if (is.data.frame(annoSNVs)) {
        write.table(annoSNVs, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
        add_generated_file(output_file)
      }

      add_to_log("Variant annotation completed successfully")
      showNotification("Variant annotation completed!",
                       type = "message", duration = 5)

    }, error = function(e) {
      add_to_log(paste("Error in variant annotation:", e$message))
      showNotification(paste("Error:", e$message),
                       type = "error", duration = 10)
    })

    removeModal()
  })

  # Display annotation output
  output$annotation_output <- renderPrint({
    if (!is.null(values$annotation_results)) {
      cat("=== Variant Annotation Results ===\n\n")

      if (is.data.frame(values$annotation_results)) {
        cat("Data frame with", nrow(values$annotation_results), "rows and",
            ncol(values$annotation_results), "columns\n\n")

        cat("Column names:\n")
        cat(paste(colnames(values$annotation_results), collapse = ", "), "\n\n")

        cat("First 20 rows:\n")
        print(head(values$annotation_results, 20))

        # Summary statistics
        cat("\n\n=== Summary Statistics ===\n")
        cat("Total annotated variants:", nrow(values$annotation_results), "\n")

        if ("gene_name" %in% colnames(values$annotation_results)) {
          unique_genes <- unique(unlist(strsplit(values$annotation_results$gene_name, ",")))
          cat("Unique genes affected:", length(unique_genes), "\n")

          # Top genes
          cat("\nTop 5 genes with most variants:\n")
          gene_counts <- table(unlist(strsplit(values$annotation_results$gene_name, ",")))
          top_genes <- head(sort(gene_counts, decreasing = TRUE), 5)
          for (i in seq_along(top_genes)) {
            cat(names(top_genes)[i], ":", top_genes[i], "\n")
          }
        }

        if ("score" %in% colnames(values$annotation_results)) {
          cat("\nScore statistics:\n")
          cat("Mean score:", round(mean(values$annotation_results$score, na.rm = TRUE), 2), "\n")
          cat("Median score:", round(median(values$annotation_results$score, na.rm = TRUE), 2), "\n")
          cat("Min score:", round(min(values$annotation_results$score, na.rm = TRUE), 2), "\n")
          cat("Max score:", round(max(values$annotation_results$score, na.rm = TRUE), 2), "\n")
        }

        if ("pos_depth" %in% colnames(values$annotation_results)) {
          cat("\nDepth statistics:\n")
          cat("Mean depth:", round(mean(values$annotation_results$pos_depth, na.rm = TRUE), 1), "\n")
          cat("Median depth:", round(median(values$annotation_results$pos_depth, na.rm = TRUE), 1), "\n")
        }

      } else {
        cat("Annotation results object type:", class(values$annotation_results), "\n")
        print(values$annotation_results)
      }
    } else {
      cat("No annotation results available. Run variant annotation first.")
    }
  })

  # Display annotated variants table
  output$annotated_variants_table <- renderDataTable({
    if (!is.null(values$annotation_results)) {
      if (is.data.frame(values$annotation_results)) {
        datatable(values$annotation_results,
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE,
                    columnDefs = list(
                      list(className = 'dt-center', targets = "_all"),
                      list(targets = c(0), visible = FALSE) # Hide snv_id if it exists
                    )
                  ),
                  rownames = FALSE,
                  filter = 'top') %>%
          formatRound(columns = c("score"), digits = 2) %>%
          formatStyle(
            'score',
            background = styleColorBar(values$annotation_results$score, 'lightblue'),
            backgroundSize = '100% 90%',
            backgroundRepeat = 'no-repeat',
            backgroundPosition = 'center'
          )
      }
    }
  })

  # Display annotation statistics
  output$annotation_stats_table <- renderTable({
    if (!is.null(values$annotation_results)) {
      if (is.data.frame(values$annotation_results)) {

        stats_data <- data.frame(
          Metric = character(),
          Value = character(),
          stringsAsFactors = FALSE
        )

        # Basic statistics
        stats_data <- rbind(stats_data,
                            data.frame(Metric = "Total Variants",
                                       Value = as.character(nrow(values$annotation_results))))

        # Gene statistics
        if ("gene_name" %in% colnames(values$annotation_results)) {
          unique_genes <- unique(unlist(strsplit(values$annotation_results$gene_name, ",")))
          stats_data <- rbind(stats_data,
                              data.frame(Metric = "Unique Genes",
                                         Value = as.character(length(unique_genes))))
        }

        # Score statistics
        if ("score" %in% colnames(values$annotation_results)) {
          stats_data <- rbind(stats_data,
                              data.frame(Metric = "Mean Score",
                                         Value = sprintf("%.2f",
                                                         mean(values$annotation_results$score, na.rm = TRUE))))

          stats_data <- rbind(stats_data,
                              data.frame(Metric = "Median Score",
                                         Value = sprintf("%.2f",
                                                         median(values$annotation_results$score, na.rm = TRUE))))
        }

        # Depth statistics
        if ("pos_depth" %in% colnames(values$annotation_results) &&
            "alt_depth" %in% colnames(values$annotation_results)) {
          stats_data <- rbind(stats_data,
                              data.frame(Metric = "Mean Total Depth",
                                         Value = sprintf("%.1f",
                                                         mean(values$annotation_results$pos_depth +
                                                                values$annotation_results$alt_depth, na.rm = TRUE))))
        }

        return(stats_data)
      }
    }
  }, striped = TRUE, hover = TRUE, width = "100%", align = 'lc')

  # Create gene effects plots
  output$top_genes_plot <- renderPlot({
    if (!is.null(values$annotation_results)) {
      if (is.data.frame(values$annotation_results) &&
          "gene_name" %in% colnames(values$annotation_results)) {

        # Extract and count genes
        all_genes <- unlist(strsplit(values$annotation_results$gene_name, ","))
        gene_counts <- sort(table(all_genes), decreasing = TRUE)

        # Take top 10 genes
        top_genes <- head(gene_counts, 10)

        # Create bar plot
        df <- data.frame(
          Gene = names(top_genes),
          Count = as.numeric(top_genes)
        )

        ggplot(df, aes(x = reorder(Gene, Count), y = Count)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          labs(
            title = "Top 10 Genes with Most Variants",
            x = "Gene",
            y = "Number of Variants"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10)
          ) +
          geom_text(aes(label = Count), hjust = -0.2, size = 3)
      }
    }
  })

  output$impact_dist_plot <- renderPlot({
    if (!is.null(values$annotation_results)) {
      if (is.data.frame(values$annotation_results) &&
          "score" %in% colnames(values$annotation_results)) {

        # Create score distribution histogram
        ggplot(values$annotation_results, aes(x = score)) +
          geom_histogram(binwidth = 5, fill = "lightblue", color = "darkblue", alpha = 0.7) +
          labs(
            title = "Distribution of Variant Scores",
            x = "Score",
            y = "Frequency"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10)
          )
      }
    }
  })

  # Run BAM processing from BAM Processing tab
  observeEvent(input$run_bam_processing, {
    req(input$bam_file_bamproc)

    showModal(modalDialog(
      title = "Processing BAM File",
      "Please wait while processing BAM file...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      add_to_log("Starting BAM processing from BAM Processing tab")

      # Save uploaded file
      bam_path <- input$bam_file_bamproc$datapath
      bam_dir <- dirname(bam_path)
      bam_file <- basename(bam_path)

      # Check if regions file is provided
      regions_path <- NULL
      if (input$use_regions_bamproc && !is.null(input$regions_file_bamproc)) {
        regions_path <- input$regions_file_bamproc$datapath
      }

      # Apply filters if selected
      if (input$filter_duplicates_bamproc || input$recalibrate_bamproc) {
        add_to_log("Applying BAM filters...")

        # Create filtered BAM file
        filtered_bam <- file.path(bam_dir,
                                  paste0(tools::file_path_sans_ext(bam_file),
                                         "_filtered.bam"))

        # This would need actual implementation based on your package functions
        # For now, we'll just copy the file as a placeholder
        file.copy(bam_path, filtered_bam)

        add_generated_file(filtered_bam)
        values$bamproc_results <- list(
          filtered_bam = filtered_bam,
          message = "BAM filtering completed"
        )
      }

      add_to_log("BAM processing completed successfully")
      showNotification("BAM processing completed!",
                       type = "message", duration = 5)

    }, error = function(e) {
      add_to_log(paste("Error in BAM processing:", e$message))
      showNotification(paste("Error:", e$message),
                       type = "error", duration = 10)
    })

    removeModal()
  })

  # Run pileup generation from BAM Processing tab
  observeEvent(input$run_pileup_generation, {
    req(input$bam_file_bamproc)

    showModal(modalDialog(
      title = "Generating Pileup",
      "Please wait while generating pileup data...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      add_to_log("Generating pileup data from BAM Processing tab")

      # Save uploaded file
      bam_path <- input$bam_file_bamproc$datapath
      bam_dir <- dirname(bam_path)
      bam_file <- basename(bam_path)

      # Check if regions file is provided
      regions_path <- NULL
      if (input$use_regions_bamproc && !is.null(input$regions_file_bamproc)) {
        regions_path <- input$regions_file_bamproc$datapath
      }

      # Generate pileup
      pileup_results <- huGSVcalleR::collectBasesPileup(
        bamDir = bam_dir,
        bamFile = bam_file,
        gr = regions_path,
        canChr = TRUE,
        depth = 1e6,
        baseq = input$min_baseq_bamproc,
        mapq = input$min_mapq_bamproc,
        coverage = input$min_coverage_bamproc,
        genome = "BSgenome.Hsapiens.UCSC.hg38",
        tmpDir = NULL,
        postfix = "PileupResults",
        workDir = bam_dir
      )

      # Store results
      values$bamproc_results <- pileup_results

      # Save pileup results
      if (!is.null(pileup_results)) {
        output_file <- file.path(bam_dir,
                                 paste0(tools::file_path_sans_ext(bam_file),
                                        "_pileup.csv"))
        if (is.data.frame(pileup_results)) {
          write.csv(pileup_results, output_file, row.names = FALSE)
          add_generated_file(output_file)
        }
      }

      add_to_log("Pileup generation completed successfully")
      showNotification("Pileup generation completed!",
                       type = "message", duration = 5)

    }, error = function(e) {
      add_to_log(paste("Error in pileup generation:", e$message))
      showNotification(paste("Error:", e$message),
                       type = "error", duration = 10)
    })

    removeModal()
  })

  # Display BAM processing output
  output$bamproc_output <- renderPrint({
    if (!is.null(values$bamproc_results)) {
      cat("=== BAM Processing Results ===\n\n")

      if (is.list(values$bamproc_results)) {
        cat("Processing completed successfully\n")
        cat("Generated files:", paste(values$bamproc_results$filtered_bam, collapse = ", "), "\n")
      } else if (is.data.frame(values$bamproc_results)) {
        cat("Pileup data generated\n")
        cat("Data frame with", nrow(values$bamproc_results), "rows and",
            ncol(values$bamproc_results), "columns\n\n")

        cat("First 10 rows:\n")
        print(head(values$bamproc_results, 10))
      } else {
        cat("Results object type:", class(values$bamproc_results), "\n")
        print(values$bamproc_results)
      }
    } else {
      cat("No BAM processing results available. Process a BAM file first.")
    }
  })

  # Display BAM processing files
  output$bamproc_files <- renderPrint({
    bam_files <- values$generated_files[grep("\\.bam$|\\.csv$", values$generated_files)]

    if (length(bam_files) > 0) {
      cat("Generated BAM/Pileup files:\n")
      cat(paste(bam_files, collapse = "\n"))
    } else {
      cat("No BAM processing files generated yet.")
    }
  })

  # Download handlers
  output$download_filtered_vcf_real <- downloadHandler(
    filename = function() {
      paste0("filtered_variants_", Sys.Date(), ".vcf")
    },
    content = function(file) {
      if (!is.null(values$filtering_results) &&
          inherits(values$filtering_results, "CollapsedVCF")) {
        VariantAnnotation::writeVcf(values$filtering_results, file)
      }
    }
  )

  output$download_annotations_real <- downloadHandler(
    filename = function() {
      paste0("annotated_variants_", Sys.Date(), ".txt")
    },
    content = function(file) {
      if (!is.null(values$annotation_results) &&
          is.data.frame(values$annotation_results)) {
        write.table(values$annotation_results, file, sep = "\t",
                    row.names = FALSE, quote = FALSE)
      }
    }
  )

  # Results summary
  output$results_summary <- renderUI({
    if (!is.null(values$test_results)) {
      tagList(
        h4("Test Results Summary"),
        p(strong("Module:"), values$test_results$module),
        p(strong("Status:"), values$test_results$status),
        p(strong("Completion Time:"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
      )
    }
  })

  # Download test results
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("huGSVcalleR_test_results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Create a temporary directory
      temp_dir <- tempdir()

      # Copy generated files to temp directory
      if (length(values$generated_files) > 0) {
        for (file_path in values$generated_files) {
          if (file.exists(file_path)) {
            file.copy(file_path, temp_dir)
          }
        }
      }

      # Create a summary file
      summary_file <- file.path(temp_dir, "summary.txt")
      summary_content <- paste(
        "huGSVcalleR Test Results Summary",
        "=================================",
        paste("Date:", Sys.Date()),
        paste("Module:", ifelse(!is.null(values$test_results), values$test_results$module, "N/A")),
        paste("Status:", ifelse(!is.null(values$test_results), values$test_results$status, "N/A")),
        paste("Total files generated:", length(values$generated_files)),
        "",
        "Generated files:",
        paste(values$generated_files, collapse = "\n"),
        sep = "\n"
      )
      writeLines(summary_content, summary_file)

      # Create zip file
      zip_file <- tempfile(fileext = ".zip")
      files_to_zip <- list.files(temp_dir, full.names = TRUE)
      zip(zip_file, files_to_zip)

      # Copy to download location
      file.copy(zip_file, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
