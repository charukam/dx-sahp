\NeedsTeXFormat{LaTeX2e}[1996/12/01]
\providecommand{\classname}{ucalgarythesis}
\ProvidesClass{\classname}[2016/02/01 v1.0]


%% Load base class using current setting for basic options.
\providecommand{\baseclass}{book}
\LoadClass{\baseclass}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Class Options
%% New option: binding (adds extra margins on binding side)
\newif\if@binding
\@bindingfalse

%% Options for binding 
\DeclareOption{binding}{\@bindingtrue}
\DeclareOption{nobinding}{\@bindingfalse}

%% New option: unbold (makes frontmatter sections unbold in TOC)
\newif\if@unbold
\@unboldfalse

%% Options for unbold
\DeclareOption{unbold}{\@bindingtrue}
\DeclareOption{nounbold}{\@bindingfalse}

%% Options for twosided or onesided printing
\DeclareOption{twoside}{\@twosidetrue}
\DeclareOption{oneside}{\@twosidefalse}

%% Options for spacing (default is 1.5 spacing)
\DeclareOption{singlespaced}{\AtEndOfClass{\singlespacingnoskip}}
\DeclareOption{onehalfspaced}{\AtEndOfClass{\onehalfspacing}}
\DeclareOption{doublespaced}{\AtEndOfClass{\doublespacing}}

%% All other options are passed to the base class directly.
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{\baseclass}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Option Execution
\ExecuteOptions{letterpaper,oneside,nobinding,nounbold,12pt}
\ExecuteOptions{doublespaced}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Process options.
\ProcessOptions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Necessary Commands
\providecommand\phantomsection{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Packages
\RequirePackage{geometry}
\RequirePackage{setspace}

\geometry{hmargin={1in,1in},vmargin={1in,1in}}  
\if@binding
    \geometry{bindingoffset=0.5in}
\else
    \geometry{bindingoffset=0in}
\fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Title Page
\newcommand{\makethesistitle}
{
    \begin{titlepage}
    \thispagestyle{empty}
    \begin{center}
    {\textbf{DESIGN AND CONSTRUCTION OF A NOVEL DIRECT EXPANSION SOLAR ASSISTED HEAT PUMP (DX-SAHP)}\vspace{5mm}\\TEAM 2 \par}
    \vfill
    \begin{flushleft}
        \begin{tabular}{ l l r }
        Project Team: \hspace{26mm} & 1. Charuka Marasinghe & \ \ \ \ \ 30023058 \\ 
                                    & 2. Dhruvi Patel        & \ \ \ \ \ 30037478 \\  
                                    & 3. Edwin Pulikkottil   & \ \ \ \ \ 30024426 \\
                                    & 4. Jessica Sambi       & \ \ \ \ \ 30016532 \\  
                                    & 5. Kerwin Virtusio     & \ \ \ \ \ 30037455 \\
                                    & 6. Nadia Elgamal       & \ \ \ \ \ 30017358 \\
                                    & & \\
        Course Instructors:         & Dr. Alex Rameriez-Serrano & \\
        Teaching Assistant:         & Danny Wong & \\
                                    & & \\
        Project Sponsors:           & 1. Dr. Aggrey Mwesigye & \\
                                    & 2. Dr. Simon Li & \\
        \end{tabular}
    \end{flushleft}
    \vfill
     Department of Mechanical and Manufacturing Engineering \\
    \vfill
    {ENME 502\\April 12, 2022\\University of Calgary \par}
    
    \begin{figure}[H]
        \centering
        \includegraphics[width=2.5cm]{images/uofc.jpg}
    \end{figure}

    \end{center}
    \end{titlepage}
    \setcounter{page}{2}
    \cleardoublepage
}

%% Commands for the thesis title page
\def\dept#1{\gdef\@dept{\uppercase{#1}}}
\def\degree#1{\gdef\@degree{\uppercase{#1}}}
\def\monthname#1{\gdef\@monthname{\uppercase{#1}}}
\def\gradyear#1{\gdef\@gradyear{#1}}
\def\thesis#1{\gdef\@thesis{\uppercase{#1}}}
\def\author#1{\gdef\@author{#1}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Front Matter Sections
\if@unbold
    \let\bforigdefault\bfdefault
    \let\oldfrontmatter\frontmatter
    \renewcommand{\frontmatter}{\oldfrontmatter\addtocontents{toc}{\let\string\bfdefault\string\mddefault}}
    \let\oldmainmatter\mainmatter
    \renewcommand{\mainmatter}{\addtocontents{toc}{\let\string\bfdefault\string\bforigdefault}\oldmainmatter}
\fi

%% Dedication Environment
\newenvironment{dedication}{\par\vspace*{.15\textheight}\centering}
{\par}

%% Redefine Table of Contents     
\renewcommand\contentsname{Table of Contents}
\let\oldtableofcontents\tableofcontents
\renewcommand{\tableofcontents}
{ 
    \begin{singlespace}
    \oldtableofcontents
    \addcontentsline{toc}{chapter}{Table of Contents}
    \phantomsection
    \end{singlespace}      
}

%% Redefine List of Figures 
\let\oldlistoffigures\listoffigures
\renewcommand{\listoffigures}
{
    \begin{singlespace}
    \oldlistoffigures
    \addcontentsline{toc}{chapter}{List of Figures}
    \phantomsection
    \end{singlespace} 
}
     
%% Redefine List of Tables
\let\oldlistoftables\listoftables
\renewcommand{\listoftables}
{
    \begin{singlespace}
    \oldlistoftables
    \addcontentsline{toc}{chapter}{List of Tables}
    \phantomsection
    \end{singlespace} 
}

%% Define addsymbol Macro for List of Symbols
\def\addsymbol #1: #2{$#1$  \> \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \parbox{5in}{#2}\\}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Bibliography
\newcommand{\thesisbibliography}
{  
    \if@unbold
        \addtocontents{toc}{\let\string\bfdefault\string\mddefault}
    \fi
    \addcontentsline{toc}{chapter}{Bibliography}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Page Styles
\if@twoside
    \def\ps@plain
    {
        \let\@mkboth\@gobbletwo
        \def\@oddhead{\reset@font\hfil\thepage}
        \let\@oddfoot\@empty
        \def\@evenhead{\reset@font\thepage\hfil}
        \let\@evenfoot\@empty
    }
\fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization
\setcounter{secnumdepth}{6}      % Depth of Sections
\thesis{thesis}                  % Default is thesis
\ps@plain                        % Default page style is plain
\onecolumn                       % Single column
\if@twoside\else\raggedbottom\fi % Ragged bottom unless twoside option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% End of cls
