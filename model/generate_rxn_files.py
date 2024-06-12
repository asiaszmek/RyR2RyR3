from lxml import etree
import argparse
import sys


def read_in_files(flist):
    roots = []
    species = set()
    specie_diff = {}
    rxn_ids = set()
    my_rxn_file = etree.Element("ReactionScheme")

    for fname_xml in flist:
        print(fname_xml)
        tree = etree.parse(fname_xml)
        roots.append(tree.getroot())
        for son in roots[-1]:
            if son.tag == "Specie":
                specie_id = son.attrib["id"]
                specie_kdiff = son.attrib["kdiff"]
                species.add(specie_id)
                if specie_id in specie_diff:
                    if specie_diff[specie_id]!= specie_kdiff:
                        print(specie_id, specie_diff[specie_id], specie_kdiff)
                else:
                    specie_diff[specie_id] = specie_kdiff
                    my_rxn_file.append(son)
            elif son.tag == "Reaction":
                if son.attrib["id"] in rxn_ids:
                    print(son.attrib["id"], " already exists", fname_xml)
                else:
                    rxn_ids.add(son.attrib["id"])
                    my_rxn_file.append(son)
    return my_rxn_file

def Parser():
    parser = argparse.ArgumentParser(description='Generate a rxn file using provided reaction files')
    parser.add_argument('input', nargs='+',
                        help='input rxn files')
    parser.add_argument('--output', default="Rxn.xml",
                        help='Where to save the reaction file')
    return parser



if __name__ == "__main__":
    fnames = []
    args = Parser().parse_args()
    for name in args.input:
        if name.endswith("xml"):
            fnames.append(name)
    if not fnames:
        sys.exit('Do specify at least one reaction file')
    output = args.output
   

    my_rxn_f = read_in_files(fnames)
    with open(output, "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    
