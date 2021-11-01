import xgboost as xgb
import os
import click
from math import ceil


def convert_model_to_c(model, ion_type, filename):
    print("Initialising model")
    bst = xgb.Booster({"nthread": 64})
    bst.load_model(model)
    bst.dump_model("{}_dump.raw.txt".format(filename))
    num_nodes = []
    mmax = 0

    print("Collecting number of trees")
    with open("{}_dump.raw.txt".format(filename)) as f:
        for row in f:
            if row.startswith("booster"):
                if row.startswith("booster[0]"):
                    mmax = 0
                else:
                    num_nodes.append(mmax + 1)
                    mmax = 0
                continue
            l = int(row.rstrip().replace(" ", "").split(":")[0])
            if l > mmax:
                mmax = l
    num_nodes.append(mmax + 1)
    forest = []
    tree = None
    b = 0
    print("Creating forest")
    with open("{}_dump.raw.txt".format(filename)) as f:
        for row in f:
            if row.startswith("booster"):
                if row.startswith("booster[0]"):
                    tree = [0] * num_nodes[b]
                    b += 1
                else:
                    forest.append(tree)
                    tree = [0] * num_nodes[b]
                    b += 1
                continue
            l = row.rstrip().replace(" ", "").split(":")
            if l[1][:4] == "leaf":
                tmp = l[1].split("=")
                tree[int(l[0])] = [-1, float(tmp[1]), -1, -1]  # !!!!
            else:
                tmp = l[1].split("yes=")
                tmp[0] = tmp[0].replace("[Features", "")
                tmp[0] = tmp[0].replace("[Feature", "")
                tmp[0] = tmp[0].replace("[f", "")
                tmp[0] = tmp[0].replace("]", "")
                tmp2 = tmp[0].split("<")
                if float(tmp2[1]) < 0:
                    tmp2[1] = 1
                tmp3 = tmp[1].split(",no=")
                tmp4 = tmp3[1].split(",")
                tree[int(l[0])] = [
                    int(tmp2[0]),
                    int(ceil(float(tmp2[1]))),
                    int(tmp3[0]),
                    int(tmp4[0]),
                ]
        forest.append(tree)

        print("Writing forest in C")
        with open("{}.c".format(filename), "w") as fout:
            fout.write(
                "float score_{}_{}(unsigned int* v){{\n".format(
                    filename, ion_type
                )
            )
            fout.write("float s = 0.;\n")
            for tt in range(len(forest)):
                fout.write(tree_to_code(forest[tt], 0, 1))
            fout.write("\nreturn s;}\n")

    os.remove("{}_dump.raw.txt".format(filename))


def tree_to_code(tree, pos, padding):
    p = "\t" * padding
    if tree[pos][0] == -1:
        if tree[pos][1] < 0:
            return p + "s = s {};\n".format(tree[pos][1])
        else:
            return p + "s = s + {};\n".format(tree[pos][1])
    return p + "if (v[{}]<{}){{\n{}}}\n{}else{{\n{}}}".format(
        tree[pos][0],
        tree[pos][1],
        tree_to_code(tree, tree[pos][2], padding + 1),
        p,
        tree_to_code(tree, tree[pos][3], padding + 1),
    )


@click.command()
@click.option("--model", help="xgboost model")
@click.option("--ion_type", help="Y or B")
@click.option("--filename", help="filename for the C model")
def main(model, ion_type, filename):
    convert_model_to_c(model, ion_type, filename)


if __name__ == "__main__":
    main()
