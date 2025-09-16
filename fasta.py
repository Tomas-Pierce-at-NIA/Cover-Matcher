# FASTA file parser that also parses out some uniprot headers
import re

class StrBuild:
    "string building helper class"
    def __init__(self):
        self.arena = []
    
    def add_str(self, s :str):
        self.arena.append(s)
    
    def build(self) -> str:
        return "".join(self.arena)

class Entry:
    "represent fasta entry"
    
    __org_code_pat = re.compile(r"OX=(\S+)")
    __gene_name_pat = re.compile(r"GN=(\S+)")
    __accession_pat = re.compile(r"\|(\S+)\|")
    __ename_pat = re.compile(r"sp\|\S+\|(\S+)\s")
    
    def __init__(self, desc :str, seq :str):
        self.desc = desc
        self.seq = seq
    
    @property
    def org_code(self):
        match = self.__org_code_pat.search(self.desc)
        return match.groups(1)[0]
    
    @property
    def gene_name(self):
        match = self.__gene_name_pat.search(self.desc)
        if match is not None:
            return match.groups(1)[0]
        else:
            return None
    
    @property
    def accession(self):
        match = self.__accession_pat.search(self.desc)
        return match.groups(1)[0]
    
    @property
    def e_name(self):
        match = self.__ename_pat.search(self.desc)
        return match.groups(1)[0]

def iter_entries(readable):
    "create generator so that we can iterate directly over FASTA entries in a file"
    builder = StrBuild()
    desc = ""
    for line in readable:
        if line.startswith('>'):
            text = builder.build()
            entry = Entry(desc, text)
            if desc != "":
                yield entry
            desc = line
            builder = StrBuild()
        else:
            linetext = line.rstrip('\n\r')
            builder.add_str(linetext)
    text = builder.build()
    entry = Entry(desc, text)
    yield entry
    


if __name__ == '__main__':
    FILENAME = "uniprotkb_human_AND_model_organism_9606_2025_08_26.fasta"
    with open(FILENAME) as handle:
        entries = iter_entries(handle)
        for entry in entries:
            print(entry.desc)
            print(entry.seq)
            print(entry.org_code)
            print(entry.gene_name)
            print(entry.accession)
            print(entry.e_name)