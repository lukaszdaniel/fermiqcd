using namespace std;

class parameter {
private:

	map<string, double> _d_rep;
	map<string, int> _i_rep;
	map<string, string> _s_rep;

	int fscanf(FILE *);

public:
	parameter(FILE *fin) {
		fscanf(fin);
	}
	;

	int fprint(FILE *);

	//  bool   defined_d(string name);
	//  bool   defined_i(string name);
	//  bool   defined_s(string name);

	parameter(char *file) {
		FILE *fin;

		if ((fin = fopen(file, "r")) == NULL) {
			cerr << "cannot open file " << string(file) << endl;
			exit(1);
		}
		fscanf(fin);
		fclose(fin);
	}

	double d(string name) {
		map<string, double>::iterator p = _d_rep.find(string(name));

		if (p == _d_rep.end()) {
			cerr << "double parameter " << name << " was not defined" << endl;
			exit(1);
		}
		return p->second;
	}
	;

	int i(string name) {
		map<string, int>::iterator p = _i_rep.find(string(name));

		if (p == _i_rep.end()) {
			cerr << "int parameter " << name << " was not defined" << endl;
			exit(1);
		}
		return p->second;
	}
	;

	string s(string name) {
		map<string, string>::iterator p = _s_rep.find(string(name));

		if (p == _s_rep.end()) {
			cerr << "string parameter " << name << " was not defined" << endl;
			exit(1);
		}
		return p->second;
	}
	;

	bool defined_d(string name) {
		map<string, double>::iterator p = _d_rep.find(string(name));
		return (p != _d_rep.end());
	}
	;

	bool defined_i(string name) {
		map<string, int>::iterator p = _i_rep.find(string(name));
		return (p != _i_rep.end());
	}
	;

	bool defined_s(string name) {
		map<string, string>::iterator p = _s_rep.find(string(name));
		return (p != _s_rep.end());
	}
	;

};

int parameter::fscanf(FILE *fin) {
	const int LENGTH = 64;
	const int LINE_LENGTH = 256;
	char line[LINE_LENGTH];

	char type[LENGTH];
	char name[LENGTH];
	char value[LENGTH];

	int line_number = 0;

	while (fgets(line, LINE_LENGTH, fin) != NULL) {
		char white_char[] = " \t\n";

		if (strspn(line, white_char) != strlen(line) and line[0] != '#') {
			line_number++;

			if (sscanf(line, "%s %s %s\n", type, name, value) == 3) {
				if (strncmp(type, "int", LENGTH) == 0)
					_i_rep[string(name)] = atoi(value);
				else if (strncmp(type, "double", LENGTH) == 0)
					_d_rep[string(name)] = atof(value);
				else if (strncmp(type, "string", LENGTH) == 0)
					_s_rep[string(name)] = string(value);
				else {
					cerr << "unknown type " << string(type) << endl;
					exit(1);
				}
			} else {
				cerr << "error reading line " << line_number << endl;
				exit(1);
			}
		}
	}
	return 0;
}

int parameter::fprint(FILE *fout) {
	for (map<string, int>::const_iterator p = _i_rep.begin(); p != _i_rep.end();
			p++) {
		fprintf(fout, "int  %s =  %d\n", (p->first).c_str(), p->second);
	}

	for (map<string, double>::const_iterator p = _d_rep.begin();
			p != _d_rep.end(); p++) {
		fprintf(fout, "double  %s =  %.12g\n", (p->first).c_str(), p->second);
	}

	for (map<string, string>::const_iterator p = _s_rep.begin();
			p != _s_rep.end(); p++) {
		fprintf(fout, "string  %s  =  %s\n", (p->first).c_str(),
				(p->second).c_str());
	}
	fflush(fout);
	return 0;
}
