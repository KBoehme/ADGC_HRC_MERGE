""" Python utils for generating latex docs """


def create_table(title, header, contents, ):
    table = """\begin{table}[H] \caption{"""+title+"""}
    \begin{center}
    \begin{tabular}{lccccc}
    \hline"""+ " & ".join(header) + """ \\ \hline""" + "\n"

    for line in contents:
        table = table + " & ".join(line) + """\\""" + "\n"
    table = table + """\end{tabular}
    \end{center}
    \label{table:full}
    \end{table}
    """
    return table

def generate_latex_doc(template. output):
