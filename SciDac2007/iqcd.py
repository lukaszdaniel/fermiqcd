#!/usr/bin/env python3
import sys
import re
import os

# Define allowed commands with their expected arguments and default values
allowed_commands = {
    "cold": [["TxXxYxZ", r"\d+x\d+x\d+x\d+", "16x4x4x4"], ["nc", r"\d+", "3"]],
    "hot": [["TxXxYxZ", r"\d+x\d+x\d+x\d+", "16x4x4x4"], ["nc", r"\d+", "3"]],
    "load": [
        ["filename", r"[a-zA-Z0-9_/.\-]*", None],
        ["precision", r"(?:float)|(?:double)", "float"],
    ],
    "lload": [
        ["filename", r"[a-zA-Z0-9_/.\-\*]*", "gauge.*.mdp"],
        ["precision", r"(?:float)|(?:double)", "float"],
    ],
    "save_partitioning": [["filename", r"[a-zA-Z0-9_/.\-]*", "partitining"]],
    "heatbath": [["beta", r"\d+(.\d+)+", "5.0"], ["steps", r"\d+", "1"]],
    "save": [["filename", r"[a-zA-Z0-9_/.\-\*]*", "gauge.*.mdp"]],
    "plaquette": [],
    "ape_smear": [
        ["alpha", r"\d(.\d+)+", "0.7"],
        ["steps", r"\d+", "20"],
        ["cooling_steps", r"\d+", "10"],
    ],
    "coulomb_gauge_fix": [
        ["boost", r"\d(.\d+)+", "1.0"],
        ["steps", r"\d+", "20"],
        ["precision", r"\d+(.\d+)+(e\-\d+)*", "1e-6"],
        ["z3", r"\d", "0"],
    ],
    "landau_gauge_fix": [
        ["boost", r"\d(.\d+)+", "1.0"],
        ["steps", r"\d+", "20"],
        ["precision", r"\d+(.\d+)+(e\-\d+)*", "1e-6"],
    ],
    "topological_charge": [
        ["filename", r"[a-zA-Z0-9_/.\-]*", "topological_charge_*.vtk"],
        ["t", r"\d", "-1"],
    ],
}

comments = {
    "cold": 'creates a cold gauge configuration, i.e. all links U(x,mu)=1\nexample: "-cold nc=3 TxXxYxZ=10x4x4x4"',
    "hot": 'creates a hot gauge configuration, i.e. all links U(x,mu)=random.SU(nc)\nexample: "-hot nc=3 TxXxYxZ=10x4x4x4"',
    "load": "loads an existing gauge configuration in mdp fermiqcd format\ndetermines lattice size from file",
    "lload": "lloads in a loop an existing gauge configuration in mdp fermiqcd format\ndetermines lattice size from file",
    "save_partitioning": "save a vtk file showing parallel partitioning information",
    "heatbath": "performs a number of heatbath steps using the WilsonGaugeAction",
    "save": "saves current gauge configuration",
    "plaquette": "compute average plaquette",
    "ape_smear": "APE smearing as defined in ref...",
    "coulomb_gauge_fix": "Coulomb Gauge Fixing algorithm as defined in ref...\n(can also do z3 fixing if z3=1 as in ref...)",
    "landau_gauge_fix": "Landau Gauge Fixing algorithm as defined in ref...",
    "topological_charge": "saves a vtk showing the topological charge as defined in ref...\n(one should ape_smear first)",
}


class Command:
    """Class to represent a command with its arguments."""

    def __init__(self, command):
        self.command = command
        self.args = {}

    def __repr__(self):
        return f"{self.command}[{repr(self.args)}]"


NC_SWITCH = """
   int nc=0;
   switch(header.bytes_per_site) {
      case 4*4*1: precision=4; nc=1; break;
      case 8*4*1: precision=8; nc=1; break;
      case 4*4*4: precision=4; nc=2; break;
      case 8*4*4: precision=8; nc=2; break;
      case 4*4*9: precision=4; nc=3; break;
      case 8*4*9: precision=8; nc=3; break;
      case 4*4*16: precision=4; nc=4; break;
      case 8*4*16: precision=8; nc=4; break;
      case 4*4*25: precision=4; nc=5; break;
      case 8*4*25: precision=8; nc=5; break;
      case 4*4*36: precision=4; nc=6; break;
      case 8*4*36: precision=8; nc=6; break;
      case 4*4*49: precision=4; nc=7; break;
      case 8*4*49: precision=8; nc=7; break;
      case 4*4*64: precision=4; nc=8; break;
      case 8*4*64: precision=8; nc=8; break;
      case 4*4*81: precision=4; nc=9; break;
      case 8*4*81: precision=8; nc=9; break;
      case 4*4*100: precision=4; nc=10; break;
      case 8*4*100: precision=8; nc=10; break;
   }
"""


def parse(argv=sys.argv[2:]):
    """Parse command-line arguments."""
    instruction = " ".join(argv)
    commands = []
    loops = []

    k = 0
    errors = []
    warnings = []

    # Parse arguments and build the command list
    for s in argv:
        if s[0] == "+":
            commands.append(Command(s[1:]))
        elif s == "{":
            loops.append(k - 1)
            commands[-1].args["begin"] = k
        elif s == "}":
            k1 = loops.pop()
            commands[k1].args["end"] = k
        else:
            try:
                a, b = s.split("=")
                commands[-1].args[a] = b
            except ValueError:
                errors.append(f"Unable to parse {s}")

        k = len(commands)

    k = 0
    # Check the first algorithm is either -cold, -hot, or -load
    if not commands or commands[0].command not in ["cold", "hot", "load"]:
        errors.append("The first algorithm must be -cold, -hot, or -load")

    # Validate each command
    for s in commands:
        if s.command == "loop":
            if "end" not in s.args:
                errors.append("+loop error: missing }")
            if "n" not in s.args:
                errors.append("+loop error: missing n=... argument")
        elif s.command not in allowed_commands:
            errors.append(f"+{s.command} error: unknown algorithm")
        else:
            d = {}
            arguments = allowed_commands[s.command]
            for a, b, c in arguments:
                d[a] = 1
                if a not in s.args:
                    if c is None:
                        errors.append(f"+{s.command} error: missing {a}=... argument")
                    else:
                        s.args[a] = c
                        warnings.append(f"+{s.command} warning: assuming default argument {a}={c}")
                elif not re.fullmatch(b, s.args[a]):
                    errors.append(f"+{s.command} error: invalid argument {a}={s.args[a]}")

            for a in s.args:
                if a not in d:
                    errors.append(f"+{s.command} error: unknown argument {a}=...")

        if s.command in ["cold", "hot"] and k > 0:
            errors.append(f"+{s.command} error: cannot be used after the first algorithm!")

        k += 1

    return instruction, commands, errors, warnings


def report(errors, warnings):
    """Print errors and warnings."""
    if errors:
        print("YOU HAVE THE FOLLOWING ERRORS:")
        for e in errors:
            print(e)
    if warnings:
        print("OK BUT SOME WARNINGS:")
        for w in warnings:
            print(w)


def generate_code(instruction, warnings, commands):
    """Generate C++ code based on the parsed commands."""
    program = "/*\n"
    program += f"    python qcd.py code {instruction}\n\n"

    # Add warnings to the program header
    for warning in warnings:
        program += f"    {warning}\n"

    program += "*/\n"
    program += '#include "fermiqcd.h"\n'
    program += "using namespace MDP;\n\n"

    program += "int main(int argc, char **argv) {\n"
    program += "   mdp.open_wormholes(argc, argv);\n"
    program += "   std::string filename;\n"
    program += "   coefficients coeff;\n"

    have_gauge = False
    indent = 1
    loops = []
    k = 0

    for s in commands:
        SPACE = "\n" + "   " * indent

        # Handle gauge field creation based on command type
        if not have_gauge:
            if s.command == "cold":
                dims = s.args.get("TxXxYxZ", "").replace("x", ",")
                program += f"   int L[] = {{{dims}}};\n"
                program += "   mdp_lattice spacetime(4, L);\n"
                program += f"   int nc = {s.args.get('nc', 3)};\n"
                program += "   gauge_field U(spacetime, nc);\n"
                program += "   set_cold(U);\n"
                have_gauge = True
            elif s.command == "hot":
                dims = s.args.get("TxXxYxZ", "").replace("x", ",")
                program += f"   int L[] = {{{dims}}};\n"
                program += "   mdp_lattice spacetime(4, L);\n"
                program += f"   int nc = {s.args.get('nc', 3)};\n"
                program += "   gauge_field U(spacetime, nc);\n"
                program += "   set_hot(U);\n"
                have_gauge = True
            elif s.command == "load":
                program += f"   mdp_field_file_header header = get_info(\"{s.args.get('filename', 'input.mdp')}\");\n"
                program += "   int L[] = {header.box[0], header.box[1], header.box[2], header.box[3]};\n"
                program += NC_SWITCH
                program += "   mdp_lattice spacetime(4, L);\n"
                program += f"   int nc = header.nc;\n"
                program += "   gauge_field U(spacetime, nc);\n"
                precision = s.args.get("precision", "double")
                if precision == "float":
                    program += f'   U.load_as_float("{s.args["filename"]}");\n'
                else:
                    program += f'   U.load_as_double("{s.args["filename"]}");\n'
                have_gauge = True
        else:
            if s.command == "loop":
                j = len(loops)
                loops.append(int(s.args.get("end", s.args.get("n", 1))))
                program += f"{SPACE}for(int i{j} = 0; i{j} < {s.args.get('n', 1)}; i{j}++) {{\n"
                indent += 1
            elif s.command == "load":
                precision = s.args.get("precision", "double")
                if precision == "float":
                    program += f"{SPACE}U.load_as_float(\"{s.args['filename']}\");\n"
                else:
                    program += f"{SPACE}U.load_as_double(\"{s.args['filename']}\");\n"
            elif s.command == "lload":
                precision = s.args.get("precision", "double")
                index = len(loops) - 1
                if precision == "float":
                    program += f"{SPACE}U.load_as_float(glob(\"{s.args['filename']}\")[i{index}]);\n"
                else:
                    program += f"{SPACE}U.load_as_double(glob(\"{s.args['filename']}\")[i{index}]);\n"
            if s.command == "save_partitioning":
                program += f"{SPACE}save_partitioning_vtk(spacetime, \"{s.args['filename']}\");\n"  # NOT IMPLEMENTED
            if s.command == "save":
                program += f"{SPACE}U.save(\"{s.args['filename']}\");\n"
            if s.command == "plaquette":
                program += (f'{SPACE}mdp << "plaquette = " << average_plaquette(U) << "\\n";\n')
            if s.command == "heatbath":
                program += f"{SPACE}coeff[\"beta\"] = {s.args['beta']};\n"
                program += f"{SPACE}WilsonGaugeAction::heatbath(U, coeff, {s.args['steps']});\n"
            if s.command == "ape_smear":
                program += f"{SPACE}ApeSmearing::smear(U, {s.args['alpha']}, {s.args['steps']}, {s.args['cooling_steps']});\n"
            if s.command == "coulomb_gauge_fix":
                program += f"{SPACE}GaugeFixing::fix(U, GaugeFixing::Coulomb, {s.args['steps']}, {s.args['precision']}, {s.args['boost']}, {s.args['z3']});\n"
            if s.command == "landau_gauge_fix":
                program += f"{SPACE}GaugeFixing::fix(U, GaugeFixing::Landau, {s.args['steps']}, {s.args['precision']}, {s.args['boost']});\n"
            if s.command == "topological_charge":
                program += f"{SPACE}{{ float tc = topological_charge_vtk(U, \"{s.args['filename']}\", {s.args['t']});\n"
                program += (f'{SPACE}  mdp << "topological_charge = " << tc << "\\n"; }}\n')
        k += 1

        # Close loops when needed
        if loops and k == loops[-1]:
            loops.pop()
            indent -= 1
            program += "\n" + "   " * indent + "}\n"

    # Final closing statements
    program += "\n\n"
    program += "   mdp.close_wormholes();\n"
    program += "   return 0;\n"
    program += "}\n"

    return program


def menu():
    """Display the main menu and handle command-line inputs."""
    if len(sys.argv) == 1:
        print("usage:")
        print('"python qcd.py help" to list all available algorithms')
        print('"python qcd.py help -heatbath" for help about the heatbath algorithm')
        print('"python qcd.py code +cold TxXxYxZ=10x4x4x4 nc=3" to create a program that makes one cold gauge configuration')
        print('"python qcd.py compile filename.cpp" to compile filename.cpp (for those who hate make like me)')
        # print '"python qcd.py submit filename.exe" to submit filename.exe to your cluster
        # print '"python qcd.py monitor pid" to monitor the executable running as pid
    elif sys.argv[1] == "code":
        instruction, commands, errors, warnings = parse()
        if errors or warnings:
            report(errors, warnings)
        if not errors:
            program = generate_code(instruction, warnings, commands)
            with open("latest.cpp", "w") as f:
                f.write(program)
            print("Saved as latest.cpp")
    elif sys.argv[1] == "compile":
        try:
            filename = sys.argv[2]
        except IndexError:
            filename = "latest.cpp"
        osx = input("Do you have an OS X 10 Mac (y/n)? ").lower()
        sse = input("Do you have a processor that supports SSE2 (y/n)? ").lower()
        pos = input("Do you have FULL POSIX support (y/n)? ").lower()
        par = input("Do you want to compile with MPI support (y/n)? ").lower()
        command = f"g++ -O2 -std=c++17 {filename} -I../Libraries -o {filename.replace('.cpp', '.exe')}"
        if osx == "y":
            command += " -DOSX"
        elif sse == "y":
            command += " -DSSE2"
        if par == "y":
            command += " -DPARALLEL"
        if pos == "n":
            command += " -DNO_POSIX"
        print(command)
        os.system(command)
    elif sys.argv[1] == "help" and len(sys.argv) < 3:
        print("Available algorithms:")
        keys = list(allowed_commands.keys())
        keys.sort()  # Sorted the keys
        for key in keys:
            print(f"   -{key}")
    elif sys.argv[1] == "help" and len(sys.argv) == 3:
        key = sys.argv[2][1:]
        if key not in allowed_commands:
            print(f"Command -{key} not supported")
        else:
            print(f"INFO FOR ALGORITHM -{key}")
            print("Description:")
            print(comments[key])
            print("Arguments:")
            for a, b, c in allowed_commands[key]:
                print(f'    "{a}"', end=" ")
                if c is None:
                    print("is required")
                else:
                    print(f"is {c} by default")
    else:
        print("I do not understand your commands")


if __name__ == "__main__":
    menu()
