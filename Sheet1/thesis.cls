%% This is my personal thesis style. It is licenced under the 
%% Creative Commons 3.0 License (http://creativecommons.org/licenses/by/3.0/)

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{thesis} [2009/12/01 My personal thesis style]

\LoadClass[fontsize=11pt,
           british,
           paper=a4,
           twoside=false, %twoside
           twocolumn=false,
           openright, %openany
           headsepline,
           footsepline,
           titlepage,
           chapterprefix=true,
           appendixprefix=true,
           parindent,
           parskip=false,
           onelinecaption,
           bigheadings,
           liststotoc,
           bibliography=openstyle,
           bibtotoc,
           idxtotoc,
           numbers=endperiod,
           listsindent,
           abstracton,
           cleardoublepage=plain,
           fleqn,
           final,
           footinclude=false]{scrbook}

\RequirePackage{pslatex}
\RequirePackage{a4wide}
\RequirePackage{fancyhdr}
\RequirePackage[Lenny]{fncychap}
\RequirePackage{ifthen}
\RequirePackage{ifpdf}
\RequirePackage[austrian,english]{babel}
\RequirePackage{graphicx}

\ifpdf
  \pdfcompresslevel=9
  \pdfpageheight=297mm
  \pdfpagewidth=210mm
  \RequirePackage[         % hyperref should be last package loaded
    pdftex,
   pdftitle={},
   pdfsubject={},
   pdfauthor={},
   pdfkeywords={},
    bookmarks,
    bookmarksnumbered,
    linktocpage,
   pagebackref,
   linkcolor=red,
   anchorcolor=red,
   citecolor=green,
   urlcolor=blue,
    pdfview={FitH},
    pdfstartview={Fit},
    pdfpagemode=UseOutlines,       % open bookmarks in Acrobat
    plainpages=false,              % avoids duplicate page number problem
    pdfpagelabels,                 % avoids duplicate page number problem
  ]{hyperref}

  \renewcommand*{\backref}[1]{}
  \renewcommand*{\backrefalt}[4]{ %}% %
  \ifcase #1 %
  (Not cited.)%
  \or
  (Cited on page~#2.)%
  \else
  (Cited on pages~#2.)%
  \fi
  }
\else
  % latex
  \DeclareGraphicsExtensions{.eps}
\fi

\setlength\parskip{\medskipamount}

\DeclareMathAlphabet{\mathpzc}{OT1}{pzc}{m}{it}
%\DeclareCaptionOption{parskip}[]{}
%\DeclareCaptionOption{parindent}[]{}

\AtBeginDocument{\pagestyle{headings}}

%\usepackage[pdftitle={},colorlinks=true]{hyperref}
\newcommand{\theHalgorithm}{\arabic{algorithm}}
\renewcommand{\thefootnote}{\fnsymbol{footnote}}

\bibliographystyle{acmtrans}

\def\prbox{\vskip-\baselineskip \hbox to \hsize{\hfill\fboxsep0pt\fbox{\phantom{\vrule width5pt height5pt depth0pt}}} \par}
\newcounter{example}[chapter]
\newenvironment{example}{\refstepcounter{example}\par\medskip\noindent%
   \textbf{Example~\thechapter.\theexample.} \rmfamily}{\nopagebreak\prbox\medskip}


\newcommand{\basedon}[1]{
  \vspace*{-10mm}
  \hfill
  \begin{minipage}{0.9\textwidth}
%    \hfill \noindent 
    \flushright 
    \small \textit{#1}
  \end{minipage}
  \vspace*{10mm}
}

%\pagestyle{headings}

\newcommand{\statdecl}{
\cleardoublepage
\vspace*{2cm}
\subsection*{Statutory Declaration}
\noindent
\textit{
I declare that I have authored this thesis independently, that I have
not used other than the declared sources / resources, and that I have
explicitly marked all material which has been quoted either literally
or by content from the used sources.
}

\vspace{1cm}

\begin{flushleft}
\begin{tabular}{lclcl}
\underline{Graz, \hspace*{4cm}} & ~\hfill & \hspace*{3cm} & ~\hfill &  \underline{\hspace*{5cm}} \\
\multicolumn{1}{c}{Place, Date}  & ~ &   & ~ &  \multicolumn{1}{c}{Signature}   \\
\end{tabular}
\end{flushleft}


\vspace{3cm}
\selectlanguage{austrian}
\subsection*{Eidesstattliche Erkl{\"a}rung}
\noindent
\textit{
Ich erkl{\"a}re an Eides statt, dass ich die vorliegende Arbeit
selbstst{\"a}ndig verfasst, andere als die angegebenen Quellen/Hilfsmittel
nicht benutzt, und die den benutzten Quellen w{\"o}rtlich und inhaltlich
entnommene Stellen als solche kenntlich gemacht habe.
}

\vspace{1cm}

\begin{flushleft}
\begin{tabular}{lclcl}
\underline{Graz, am\hspace*{3.6cm}} & ~\hfill & \hspace*{3cm} & ~\hfill &  \underline{\hspace*{5cm}} \\
\multicolumn{1}{c}{Ort, Datum}  & ~ &   & ~ &  \multicolumn{1}{c}{Unterschrift}   \\
\end{tabular}
\end{flushleft}
\selectlanguage{english}
}

\newcommand*{\@supervisor}{}
\newcommand{\supervisor}[1]{\gdef\@supervisor{#1}}
\newcommand*{\@secondreviewer}{}
\newcommand{\secondreviewer}[1]{\gdef\@secondreviewer{#1}}
\newcommand*{\@institute}{}
\newcommand{\institute}[1]{\gdef\@institute{#1}}
\newcommand*{\@logo}{}
\newcommand{\logo}[1]{\gdef\@logo{#1}}


\renewcommand*{\maketitle}{%
  \cleardoublepage
  \thispagestyle{empty}
  \null
  \vskip 2em%
  \begin{center}%
    \ifx\@subject\@empty \else
      {\subject@font \@subject \par}
      \vskip 5em
    \fi
    {\hrule\vskip 2em \titlefont\huge\@title\vspace*{1em}\hrule\vspace*{2em}\par}%
    \vskip 1.5em%
    {\Large
      \lineskip .5em%
      \begin{tabular}[t]{c}
        \@author
      \end{tabular}\par
    }%
    {\Large \@date \par}%
    \vskip \z@ \@plus 1em%
    \vfill
    {\normalsize \it \@institute \par \vspace{-2mm}}%
    {\normalsize \it Graz University of Technology \par}%
    \vskip 1em%
    {\includegraphics[width=3cm,keepaspectratio]{\@logo} \par}%
    \vspace{20mm}%
    \ifx\@supervisor\@empty 
      \vskip 2em 
    \else
      {\normalsize Supervisor/First reviewer: \@supervisor \par}%
    \fi
    \vspace{-2mm}%
    \ifx\@secondreviewer\@empty 
      \vskip 2em 
    \else
      {\normalsize Second reviewer: \@secondreviewer \par}%
    \fi
  \end{center}%
  \par
  \vskip 2em
  \next@tpage
  \null
}


