%%
%% This is file `csiroreport2017.cls'.
%%
%% IMPORTANT NOTICE:
%% If you can improve this file then please do so and send revisions to brand@csiro.au
%%
%%
\NeedsTeXFormat{LaTeX2e}[1994/12/01] \ProvidesClass{csiroreport2017}
              [2017/02/07, v1.4]

\newif\if@times
\newif\if@mtpro
\newif\if@colortbl
\newif\if@fleqn

\newif\if@review
\newif\if@final

\newif\if@formal
\newif\if@standard

\@formalfalse \DeclareOption{formal}{\@formaltrue}
\@standardfalse \DeclareOption{standard}{\@standardtrue}

\newif\if@blueberry
\newif\if@forest
\newif\if@midday
\newif\if@mint
\newif\if@sky
\newif\if@teal

\@blueberryfalse \DeclareOption{blueberry}{\@blueberrytrue}
\@forestfalse \DeclareOption{forest}{\@foresttrue}
\@middayfalse \DeclareOption{midday}{\@middaytrue}
\@mintfalse \DeclareOption{mint}{\@minttrue}
\@skyfalse \DeclareOption{sky}{\@skytrue}
\@tealfalse \DeclareOption{teal}{\@tealtrue}

\@reviewfalse \DeclareOption{review}{\@reviewtrue}

\@finalfalse \DeclareOption{final}{\@finaltrue}

\@timesfalse \DeclareOption{times}{\@timestrue}

\@mtprofalse \DeclareOption{mtpro}{\@mtprotrue}

\@colortblfalse \DeclareOption{colortbl}{\@colortbltrue}

\@fleqnfalse \DeclareOption{fleqn}{\@fleqntrue}


\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}} \ProcessOptions
\LoadClass{article}

\if@times \RequirePackage[scaled=0.90]{helvet} \RequirePackage{times} \fi

\if@fleqn
  \RequirePackage{fleqn}
  \setlength{\mathindent}{4em}
\fi

 \RequirePackage{wallpaper}

 \textwidth=160mm
 \textheight=240mm
 \topsep=0ex
 \headheight=0ex
 \headsep=0ex
 \topmargin=0mm
 \footskip=10mm
 \oddsidemargin=6pt
 \evensidemargin=0pt\parskip=1.0ex
 \parindent=0em
 \raggedbottom


\marginparwidth=35mm
\marginparsep=1ex

\RequirePackage{ifthen} \RequirePackage{parskip}


\if@mtpro
  \RequirePackage{mtpro}
\fi

\newcounter{pdfout}


%% latex or pdflatex?
\@ifundefined{pdfoutput}{  %latex-mode
  \setcounter{pdfout}{0}
  \RequirePackage{graphicx}
  \RequirePackage{color}
  \RequirePackage[%draft,
    dvips,
    pdftitle={CSIRO LaTeX report class},
    pdfauthor={brand@csiro.au},
    pdfsubject={Copyright CSIRO},
    pdfcreator={LaTeX source to pdf},
    bookmarksnumbered,bookmarksopen,pdfstartview=Fit,breaklinks=true,
    %colorlinks=true,urlcolor=blue, citecolor = magenta, linkcolor=red]{hyperref}
    colorlinks=true,urlcolor=dark-colour, citecolor = black, linkcolor=black]{hyperref}
}
{  %pdflatex-mode
  \setcounter{pdfout}{1}
  \RequirePackage{graphicx}
  \RequirePackage{color}
  \RequirePackage[
    pdftitle={CSIRO LaTeX report class},
    pdfauthor={brand@csiro.au},
    pdfsubject={Copyright CSIRO},
    pdfcreator={LaTeX source to pdf},
    bookmarksnumbered,bookmarksopen,pdfstartview=Fit,
    %colorlinks=true,urlcolor=blue, citecolor = magenta, linkcolor=red]{hyperref}
    colorlinks=true,urlcolor=dark-colour, citecolor = black, linkcolor=black]{hyperref}
}

\if@colortbl \RequirePackage{colortbl} \fi

\newboolean{formal}
\setboolean{formal}{false}
\newboolean{standard}
\setboolean{standard}{false}

\if@formal
\setboolean{formal}{true}
\fi
\if@standard
\setboolean{standard}{true}
\fi


\newboolean{blueberry}
\setboolean{blueberry}{false}
\newboolean{forest}
\setboolean{forest}{false}
\newboolean{midday}
\setboolean{midday}{false}
\newboolean{mint}
\setboolean{mint}{false}
\newboolean{sky}
\setboolean{sky}{false}
\newboolean{teal}
\setboolean{teal}{false}

\if@blueberry
\setboolean{blueberry}{true}
\fi
\if@forest
\setboolean{forest}{true}
\fi
\if@midday
\setboolean{midday}{true}
\fi
\if@mint
\setboolean{mint}{true}
\fi
\if@sky
\setboolean{sky}{true}
\fi
\if@teal
\setboolean{teal}{true}
\fi


\ifthenelse{\boolean{blueberry}}{
 \definecolor{light-colour}{rgb}{0.624,0.682,0.898} %R159 G174 B229 CSIRO Lavender
 \definecolor{dark-colour}{rgb}{0.118,0.133,0.667} %R30 G34 B170 CSIRO Blueberry
 }
 {
 \ifthenelse{\boolean{forest}}{
  \definecolor{light-colour}{rgb}{0.471,0.745,0.125} %R120 G190 B32 CSIRO Light Forest
  \definecolor{dark-colour}{rgb}{0.267,0.412,0.239} %R68 G105 B61 CSIRO Dark Forest
   }
   {
   \ifthenelse{\boolean{midday}}{
    \definecolor{light-colour}{rgb}{0,0.663,0.808} %R0 G169 B206 CSIRO Midday Blue
    \definecolor{dark-colour}{rgb}{0,0.192,0.235} %R0 G49 B60 CSIRO Midnight Blue
    }
    {
    \ifthenelse{\boolean{mint}}{
     \definecolor{light-colour}{rgb}{0.443,0.800,0.596} %R113 G204 B152 CSIRO Light Mint
     \definecolor{dark-colour}{rgb}{0,0.478,0.325} %R0 G122 B83 CSIRO Dark Mint
     }
     {
     \ifthenelse{\boolean{sky}}{
      \definecolor{light-colour}{rgb}{0.255,0.714,0.902} %R65 G182 B230 CSIRO Sky Blue
      \definecolor{dark-colour}{rgb}{0,0.294,0.529} %R0 G75 B135 CSIRO Ocean Blue
      }
      {
      \ifthenelse{\boolean{teal}}{
       \definecolor{light-colour}{rgb}{0.176,0.800,0.827} %R45 G204 B211 CSIRO Light Teal
       \definecolor{dark-colour}{rgb}{0,0.451,0.467} %R0 G115 B119 CSIRO Dark Teal
       }
       {
        % Make midday the default
        \definecolor{light-colour}{rgb}{0,0.663,0.808} %R0 G169 B206 CSIRO Midday Blue
        \definecolor{dark-colour}{rgb}{0,0.192,0.235} %R0 G49 B60 CSIRO Midnight Blue
       }
      }
     }
    }
  }
 }


% The following commands redefine the styles for the different section headings.
\makeatletter
\renewcommand{\section}{\resetequation\resettable\resetfigure%
  \@startsection{section}{1}{0em}%
  {-\baselineskip}{0.5ex}{\normalfont\Large\bfseries\sffamily\color{dark-colour}}}
\renewcommand{\subsection}{\@startsection{subsection}{2}{0em}%
  {-0.5\baselineskip}{0.1ex}{\normalfont\large\bfseries\sffamily\color{dark-colour}}}
\renewcommand{\subsubsection}{\@startsection{subsubsection}{3}{0em}%
  {-0.1\baselineskip}{0.1ex}{\normalfont\bfseries\sffamily\color{dark-colour}}}
\renewcommand{\paragraph}{\@startsection{paragraph}{4}{0em}%
  {-0.1\baselineskip}{0.01ex}{\normalfont\small\bfseries\sffamily\color{dark-colour}}}
\renewcommand{\subparagraph}{\@startsection{subparagraph}{5}{0em}%
  {-0.1\baselineskip}{0.01ex}{\normalfont\small\sffamily}\color{dark-colour}}
\renewcommand{\l@section}{\bigskip\bfseries\@dottedtocline{1}{1.0em}{1.2em}}
\renewcommand{\l@subsection}{\smallskip\normalfont\@dottedtocline{2}{2.2em}{2.2em}}
\renewcommand{\l@subsubsection}{\normalfont\@dottedtocline{3}{4.4em}{2.9em}}
\renewcommand{\l@paragraph}{\normalfont\@dottedtocline{3}{7.3em}{3.3em}}

\renewcommand{\contentsline}[4]{\sffamily\color{dark-colour}\csname l@#1\endcsname{\sffamily #2}{\color{dark-colour}\sffamily#3}{}}


\makeatother

\renewcommand{\theequation}{\thesection.\arabic{equation}}
\newcommand{\resetequation}{\setcounter{equation}{0}}
\renewcommand{\thetable}{\thesection.\arabic{table}}
\newcommand{\resettable}{\setcounter{table}{0}}
\renewcommand{\thefigure}{\thesection.\arabic{figure}}
\newcommand{\resetfigure}{\setcounter{figure}{0}}


\setcounter{secnumdepth}{4} \setcounter{tocdepth}{4}

  \newcommand{\thedoctitle}{}
  \newcommand{\doctitle}[1][]{\renewcommand{\thedoctitle}{#1}}

  \newcommand{\thedocsubtitle}{}
  \newcommand{\docsubtitle}[1][]{\renewcommand{\thedocsubtitle}{#1}}


  \newcommand{\thecommercial}{}
  \newcommand{\docinconfidence}[1][]{\renewcommand{\thecommercial}{#1}}

  \newcommand{\theabstract}{}
  \newcommand{\docabstract}[1][]{\renewcommand{\theabstract}{#1}}

  \newcommand{\thedivision}{}
  \newcommand{\docdivision}[1][]{\renewcommand{\thedivision}{#1}}

  \newcommand{\theauthors}{}
  \newcommand{\docauthors}[1][]{\renewcommand{\theauthors}{#1}}

  \newcommand{\thereportnum}{}
  \newcommand{\docreportnum}[1][]{\renewcommand{\thereportnum}{#1}}

  \newcommand{\thenumclientcopy}{0}
  \newcommand{\numclientcopy}[1][]{\renewcommand{\thenumclientcopy}{#1}}

  \newcommand{\thenumlibcopy}{0}
  \newcommand{\numlibcopy}[1][]{\renewcommand{\thenumlibcopy}{#1}}

  \newcommand{\thenumslcopy}{0}
  \newcommand{\numslcopy}[1][]{\renewcommand{\thenumslcopy}{#1}}

  \newcommand{\thenumauthorcopy}{0}
  \newcommand{\numauthorcopy}[1][]{\renewcommand{\thenumauthorcopy}{#1}}

  \newcommand{\thereportdate}{}
  \newcommand{\docreportdate}[1][]{\renewcommand{\thereportdate}{#1}}

  \newcommand{\thecopyrightyear}{}
  \newcommand{\doccopyrightyear}[1][]{\renewcommand{\thecopyrightyear}{#1}}

  \newcommand{\thefootertitle}{}
  \newcommand{\docfootertitle}[1][]{\renewcommand{\thefootertitle}{#1}}

  \newcommand{\thereportday}{}
  \newcommand{\docreportday}[1][]{\renewcommand{\thereportday}{#1\ }{}}

  \newcommand{\thereportmonth}{}
  \newcommand{\docreportmonth}[1][]{\renewcommand{\thereportmonth}{#1\ }}

  \newcommand{\thereportyear}{}
  \newcommand{\docreportyear}[1][]{\renewcommand{\thereportyear}{#1}}

  \newcommand{\theclient}{}
  \newcommand{\docclient}[1][]{\renewcommand{\theclient}{#1}}

  \newcommand{\theclientcontact}{}
  \newcommand{\docclientcontact}[1][]{\renewcommand{\theclientcontact}{#1}}

  \newcommand{\thebusinessunit}{}
  \newcommand{\docbusinessunit}[1][]{\renewcommand{\thebusinessunit}{#1}}


  \newcommand{\thefurtherinfoA}{}
  \newcommand{\docfurtherinfoA}[6][]{\renewcommand{\thefurtherinfoA}{%

  FOR FURTHER INFORMATION

   \vspace{1ex}

  \textbf{#1}\\
  #2\\
  \begin{tabular}{ll}
  \textbf{t} & #3 \\
  \textbf{e} &  \href{mailto: #4}{\ifthenelse{\boolean{formal}}{#4}
  {\ifthenelse{\boolean{standard}}{\color{white}#4}{#4}}}\\
  \textbf{w} &  \href{#5}{\ifthenelse{\boolean{formal}}{#6}
  {\ifthenelse{\boolean{standard}}{\color{white}#6}{#6}}}\\
  \end{tabular}}}

  \newcommand{\thefurtherinfoB}{}
  \newcommand{\docfurtherinfoB}[6][]{\renewcommand{\thefurtherinfoB}{%
  \textbf{#1}\\
  #2\\
  \begin{tabular}{ll}
  \textbf{t} & #3 \\
  \textbf{e} &  \href{mailto: #4}{\ifthenelse{\boolean{formal}}{#4}
  {\ifthenelse{\boolean{standard}}{\color{white}#4}{#4}}}\\
  \textbf{w} &  \href{#5}{\ifthenelse{\boolean{formal}}{#6}
  {\ifthenelse{\boolean{standard}}{\color{white}#6}{#6}}}\\
  \end{tabular}}}

  \newcommand{\thefurtherinfoC}{}
  \newcommand{\docfurtherinfoC}[6][]{\renewcommand{\thefurtherinfoC}{%
  \textbf{#1}\\
  #2\\
  \begin{tabular}{ll}
  \textbf{t} & #3 \\
  \textbf{e} &  \href{mailto: #4}{\ifthenelse{\boolean{formal}}{#4}
  {\ifthenelse{\boolean{standard}}{\color{white}#4}{#4}}}\\
  \textbf{w} &  \href{#5}{\ifthenelse{\boolean{formal}}{#6}
  {\ifthenelse{\boolean{standard}}{\color{white}#6}{#6}}}\\
  \end{tabular}}}




\newcommand{\thefrontpicture}{}
\newcommand{\docfrontpicture}[2][scale=1.0]{\renewcommand{\thefrontpicture}{%
\AddToShipoutPicture*{%
\AtPageLowerLeft{%
\parbox[b][\paperheight]{\paperwidth}{%
\hspace{-28.5mm}\raisebox{150.7mm}[0mm][0mm]{%
\begin{minipage}[c]{300mm}%
\centering%
\includegraphics[#1]{#2}%
\end{minipage}}%
}}}}}


\newcounter{partcount}
\setcounter{partcount}{1}

\newcommand{\partpage}[1]{
\clearpage

\thispagestyle{empty}

\ifthenelse{\boolean{blueberry}}{
 \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportpart-blueberry}}
 {\ifthenelse{\boolean{forest}}{
   \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportpart-forest}}
   {\ifthenelse{\boolean{midday}}{
    \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportpart-midday}}
    {\ifthenelse{\boolean{mint}}{
     \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportpart-mint}}
     {\ifthenelse{\boolean{sky}}{
      \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportpart-sky}}
      {\ifthenelse{\boolean{teal}}{
       \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportpart-teal}}
       {% Make midday the default
        \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportpart-midday}}
      }
     }
    }
   }
  }


\vspace{20mm}

{\Huge\color{light-colour}\sffamily\bfseries Part \Roman{partcount}\quad \color{dark-colour} #1}


\addcontentsline{toc}{part}{Part \Roman{partcount}\quad
#1}

\refstepcounter{partcount}

\clearpage

}


\newcommand{\thelogoA}{}
\newcommand{\doclogoA}[4]{\renewcommand{\thelogoA}{%
\AddToShipoutPicture*{%
\AtPageLowerLeft{%
\parbox[b][\paperheight]{\paperwidth}{%
\raisebox{#2}[0mm][0mm]{%
\hspace{#1}\includegraphics[#3]{#4}%
} }}}}}

\newcommand{\thelogoB}{}
\newcommand{\doclogoB}[4]{\renewcommand{\thelogoB}{%
\AddToShipoutPicture*{%
\AtPageLowerLeft{%
\parbox[b][\paperheight]{\paperwidth}{%
\raisebox{#2}[0mm][0mm]{%
\hspace{#1}\includegraphics[#3]{#4}%
} }}}}}


\newcommand{\thelogoC}{}
\newcommand{\doclogoC}[4]{\renewcommand{\thelogoC}{%
\AddToShipoutPicture*{%
\AtPageLowerLeft{%
\parbox[b][\paperheight]{\paperwidth}{%
\raisebox{#2}[0mm][0mm]{%
\hspace{#1}\includegraphics[#3]{#4}%
} }}}}}


\newcommand{\setreportfront}[1]{
\ifthenelse{\boolean{blueberry}}{
 \ClearWallPaper\thefrontpicture\ThisCenterWallPaper{1.005}{figures/reportfront-#1-blueberry}}
 {\ifthenelse{\boolean{forest}}{
   \ClearWallPaper\thefrontpicture\ThisCenterWallPaper{1.005}{figures/reportfront-#1-forest}}
   {\ifthenelse{\boolean{midday}}{
    \ClearWallPaper\thefrontpicture\ThisCenterWallPaper{1.005}{figures/reportfront-#1-midday}}
    {\ifthenelse{\boolean{mint}}{
     \ClearWallPaper\thefrontpicture\ThisCenterWallPaper{1.005}{figures/reportfront-#1-mint}}
     {\ifthenelse{\boolean{sky}}{
      \ClearWallPaper\thefrontpicture\ThisCenterWallPaper{1.005}{figures/reportfront-#1-sky}}
      {\ifthenelse{\boolean{teal}}{
       \ClearWallPaper\thefrontpicture\ThisCenterWallPaper{1.005}{figures/reportfront-#1-teal}}
       {% Make midday the default
        \ClearWallPaper\thefrontpicture\ThisCenterWallPaper{1.005}{figures/reportfront-#1-midday}}
      }
     }
    }
   }
  }
}




\newcommand{\setfronttext}{
 \AtPageLowerLeft{
  \parbox[b][\paperheight]{\paperwidth}{
   \vspace{19mm}%
   \hspace{21mm}\begin{minipage}[t]{135mm}
   \raggedright\color{white}

   \sffamily\bfseries

   \thedivision

   \ifthenelse{\boolean{standard}}{
   \color{light-colour}}{\color{dark-colour}}

   www.csiro.au

   \end{minipage} \vfill
  }
 }

 \AtPageLowerLeft{
  \parbox[b][\paperheight]{\paperwidth}{
   \vspace{56mm}%
   \hspace{21mm}\begin{minipage}[t]{180mm}
   \raggedright


   \ifthenelse{\boolean{standard}}{
   \color{white}}{\color{light-colour}}

   \Large\sffamily\bfseries\thedoctitle

   \ifthenelse{\boolean{standard}}{
   \color{white}}{\color{dark-colour}}

   \large\sffamily\bfseries\thedocsubtitle


   \ifthenelse{\boolean{standard}}{
   \color{white}}{\color{black}}

   \large\theauthors

   \large\normalfont\sffamily\thereportnum

   \thereportdate

  \vspace{10mm}

   \theclient

   \theclientcontact

   \end{minipage} \vfill
  }
 }
 \AtPageLowerLeft{
  \parbox[b][\paperheight]{\paperwidth}{
   \vspace{170mm}
   \hspace{21mm}\begin{minipage}[t]{140mm}
   \raggedright

   \ifthenelse{\boolean{standard}}{
   \color{white}}{\color{black}}

   \sffamily\large\thecommercial
   \end{minipage}\vfill
  }
 }

}

\newcommand{\setreportback}[1]{
\ifthenelse{\boolean{blueberry}}{
 \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportback-#1-blueberry}}
 {\ifthenelse{\boolean{forest}}{
   \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportback-#1-forest}}
   {\ifthenelse{\boolean{midday}}{
    \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportback-#1-midday}}
    {\ifthenelse{\boolean{mint}}{
     \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportback-#1-mint}}
     {\ifthenelse{\boolean{sky}}{
      \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportback-#1-sky}}
      {\ifthenelse{\boolean{teal}}{
       \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportback-#1-teal}}
       {% Make midday the default
        \mbox{}\ClearWallPaper\ThisCenterWallPaper{1.005}{figures/reportback-#1-midday}}
      }
     }
    }
   }
  }
}

\newcommand{\contactcsiro}{
 \AtPageLowerLeft{
  \parbox[b][\paperheight]{\paperwidth}{
   \vspace{55mm}%
   \hspace{21mm}\begin{minipage}[t]{60mm}
   \raggedright

   \small\sffamily

   \ifthenelse{\boolean{standard}}{
   \color{white}}{}

   CONTACT US

   \vspace{1ex}

   \begin{tabular}{@{}ll}
   \textbf{t} & 1300 363 400 \\

    & +61 3 9545 2176\\

   \textbf{e} &  \href{mailto: csiroenquiries@csiro.au}{
   \ifthenelse{\boolean{standard}}{{\color{white}csiroenquiries@csiro.au}}{csiroenquiries@csiro.au}}\\
   \textbf{w} &  \href{http://www.csiro.au}{
   \ifthenelse{\boolean{standard}}{{\color{white}www.csiro.au}}{www.csiro.au}}\\
   \end{tabular}

   \vspace{2ex}

   WE DO THE EXTRAORDINARY EVERY DAY 

   We innovate for tomorrow and help improve today – for our customers, 
   all Australians and the world. 

   Our innovations contribute billions of dollars to the Australian economy 
   every year. As the largest patent holder in the nation, our vast wealth 
   of intellectual property has led to more than 150 spin-off companies. 

   With more than 5,000 experts and a burning desire to get things done, 
   we are Australia’s catalyst for innovation. 

   WE IMAGINE. WE COLLABORATE. WE INNOVATE.


   \end{minipage} \vfill
  }
 }
}

\newcommand{\furtherinfo}{
 \AtPageLowerLeft{
  \parbox[b][\paperheight]{\paperwidth}{
   \vspace{55mm}%
   \hspace{90mm}\begin{minipage}[t]{80mm}
   \raggedright
   \small\sffamily
   \ifthenelse{\boolean{standard}}{\color{white}}{}

   \thefurtherinfoA

   \vspace{10mm}

   \thefurtherinfoB

   \vspace{10mm}

   \thefurtherinfoC

   \end{minipage} \vfill
  }
 }
}





\makeatletter
\renewcommand{\@evenfoot}{{\thepage~$|$~\thefootertitle}\hfil}
\renewcommand{\@oddfoot}{\hfil{\thefootertitle~$|$~\thepage}}
\makeatother


\renewcommand{\familydefault}{\sfdefault}


%=================================================

\AtBeginDocument{

  \pagenumbering{alph}

  \thispagestyle{empty}


  \ifthenelse{\boolean{formal}}{
    \setreportfront{formal}

    \mbox{}

    \AddToShipoutPicture*{
      \setfronttext
    }
  }
  {\ifthenelse{\boolean{standard}}{
    \setreportfront{standard}

    \mbox{}

    \AddToShipoutPicture*{
      \setfronttext
    }
  }
  {% Make the default front the same as formal
    \setreportfront{formal}

    \mbox{}

    \AddToShipoutPicture*{
      \setfronttext
    }
  }
  }

  \renewcommand{\thefrontpicture}{}
  \thelogoA  \renewcommand{\thelogoA}{}
  \thelogoB  \renewcommand{\thelogoB}{}
  \thelogoC  \renewcommand{\thelogoC}{}

  \clearpage

  \thispagestyle{empty}

  \thebusinessunit

  \mbox{}

  \vspace{10mm}

  {\sffamily\bfseries\color{dark-colour} Copyright and disclaimer}

  \smallskip

  \copyright\ \thecopyrightyear\ CSIRO To the extent permitted by law, all
  rights are reserved and no part of this publication covered by copyright may
  be reproduced or copied in any form or by any means except with the written
  permission of CSIRO.

  \vspace{10mm}

  {\sffamily\bfseries\color{dark-colour} Important disclaimer}

  \smallskip

  CSIRO advises that the information contained in this publication comprises
  general statements based on scientific research. The reader is advised and
  needs to be aware that such information may be incomplete or unable to be
  used in any specific situation. No reliance or actions must therefore be made
  on that information without seeking prior expert professional, scientific and
  technical advice. To the extent permitted by law, CSIRO (including its
  employees and consultants) excludes all liability to any person for any
  consequences, including but not limited to all losses, damages, costs,
  expenses and any other compensation, arising directly or indirectly from
  using this publication (in part or in whole) and any information or material
  contained in it.

  \clearpage

  \pagenumbering{roman}

  \setcounter{page}{1}

  % set up the table of contents, list of figures, tables, etc if needed.
   \parskip=0ex
   \tableofcontents

   \clearpage

   \ifthenelse{\isodd{\value{page}}}
  {}
  {\mbox{}\clearpage}


   \parskip=1ex
   \parindent=0em

  \pagenumbering{arabic}
  \setcounter{page}{1}

}

%=================================================

\AtEndDocument{

% This sets up the back cover. Whether you need the extra \clearpage will depend
% on whether the document is to be printed as double sided and whether the last page
% of text is an even or odd numbered page.

  \clearpage

  \ifthenelse{\isodd{\value{page}}}%
  {\pagestyle{empty}\mbox{}\clearpage\mbox{}}%
  {\pagestyle{empty}\mbox{}\clearpage\mbox{}\clearpage}

  \ifthenelse{\boolean{formal}}{
    \setreportback{formal}
    \AddToShipoutPicture*{
      \contactcsiro
    }
    \AddToShipoutPicture*{
      \furtherinfo
    }
  }
  {
    \ifthenelse{\boolean{standard}}{
      \setreportback{standard}
      \AddToShipoutPicture*{
        \contactcsiro
      }
      \AddToShipoutPicture*{
        \furtherinfo
      }
    }
  {
    % Make the default back the same as formal
      \setreportback{formal}
      \AddToShipoutPicture*{
        \contactcsiro
      }
      \AddToShipoutPicture*{
        \furtherinfo
      }
    }
  }

}

%=================================================

\endinput
