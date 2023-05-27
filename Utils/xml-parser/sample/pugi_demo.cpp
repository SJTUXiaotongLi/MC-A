#include <iostream>
#include "../pugixml/pugixml.hpp"

using namespace std;
using namespace pugi;

int main()
{
	cout << "Parsing a parameter file (sample.xml).....\n";
	
	xml_document doc;
	int L0, L1, N_TemperaturePoints;
	
	// load the XML file
	if ( ! doc.load_file("sample.xml") ) {
		cerr << "sample.xml is not found" << endl;
		return -1;
	}
	
	if (1) {
		xml_node tools = doc.child("ParameterData").child("System");
		for (pugi::xml_node tool: tools.children("Parameter")) {
			std::cout << "Parameter: ";
			for (pugi::xml_attribute attr: tool.attributes()) {
				std::cout << " " << attr.name() << " ... " << attr.value() << std::endl;
				if (attr.name() == std::string("L0")) {
					L0 = attr.as_int();
					cout << "L0 = " << L0 << endl;
				}
				if (attr.name() == std::string("L1")) {
					L1 = attr.as_int();
					cout << "L1 = " << L1 << endl;
				}
			}
		}
	}
	
	if (1) {
		xml_node tools = doc.child("ParameterData").child("Model");
		for (pugi::xml_node tool: tools.children("Parameter")) {
			std::cout << "Parameter: ";
			for (pugi::xml_attribute attr: tool.attributes()) {
				std::cout << " " << attr.name() << " ... " << attr.value() << std::endl;
				if (attr.name() == std::string("N_TemperaturePoints")) {
					N_TemperaturePoints = attr.as_int();
					cout << "N_TemperaturePoints = " << N_TemperaturePoints << endl;
				}
			}
		}
	}
	return 0;
}
