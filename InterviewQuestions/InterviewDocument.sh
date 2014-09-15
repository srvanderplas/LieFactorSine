% { shopt -s nullglob;
  for file in ./*.png;
    do echo "\begin{minipage}[b]{.5\linewidth}\begin{center}\includegraphics[keepaspectratio=TRUE]{$file}\captionof{figure}{\underline{\hspace{2in}}}\end{center}\end{minipage}";
  done;
}
