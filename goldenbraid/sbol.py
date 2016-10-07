import re
import xml.etree.ElementTree as ET
from xml.dom import minidom

from goldenbraid.models import Feature

HOMESPACE = "https://gbcloning.upv.es/accession"

GOLDEN_TO_SO = {}


def _parse_parts(seqrecord):
    if seqrecord.description == '<unknown description>':
        return []
    return [p for p in re.split('\(|\)|,', seqrecord.description) if p]


def convert_to_sbol(seqrecord):
    root = ET.Element('rdf:RDF')
    root.set('xmlns:dcterms', "http://purl.org/dc/terms/#")
    root.set('xmlns:prov', "http://www.w3.org/ns/prov#")
    root.set('xmlns:rdf', "http://www.w3.org/1999/02/22-rdf-syntax-ns#")
    root.set('xmlns:sbol', "http://sbols.org/v2#")
    # Main_comp
    sequence_xml = ET.SubElement(root, 'sbol:Sequence')
    sequence_xml = populate_sequence_element(sequence_xml, seqrecord)
    components = []
    for part in _parse_parts(seqrecord):
        comp_def = ET.SubElement(root, 'sbol:ComponentDefinition')
        populate_comp_def(comp_def, name=part)
        components.append(comp_def.attrib['rdf:about'])

    main_part = ET.SubElement(root, 'sbol:ComponentDefinition')
    populate_comp_def(main_part, name=seqrecord.id, components=components,
                      sequence_xml=sequence_xml)

    root_string = ET.tostring(root, 'utf-8')
    reparsed = minidom.parseString(root_string)
    return reparsed.toprettyxml(indent="  ")


def populate_comp_def(comp_def_xml, name, components=None, sequence_xml=None):
    name = name
    uri = HOMESPACE + '/' + name
    comp_def_xml.set('rdf:about', uri)
    displayid = ET.SubElement(comp_def_xml, 'sbol:displayId')
    displayid.text = name
    persisten_ident = ET.SubElement(comp_def_xml, 'sbol:persistentIdentity')
    persisten_ident.set('rdf:resource', uri)
    ET.SubElement(comp_def_xml, 'sbol:type', {'rdf:resource':
                  'http://www.biopax.org/release/biopax-level3.owl#DnaRegion'})
    if sequence_xml is not None:
        sequence = ET.SubElement(comp_def_xml, 'sbol:sequence')
        seq_xml_uri = sequence_xml.attrib['rdf:about']
        sequence.set('rdf:resource', seq_xml_uri)

    compliant_component_uris = []
    if components:
        components_xml = ET.SubElement(comp_def_xml, 'sbol:component')
        for component_uri in components:
            component_displayid = component_uri.split('/')[-1]
            compliant_component_uri = uri + '/' + component_displayid
            compliant_component_uris.append(compliant_component_uri)
            comp_xml = ET.SubElement(components_xml, 'sbol:Component')
            comp_xml.set('rdf:about', compliant_component_uri)
            sbol_displayid = ET.SubElement(comp_xml, 'sbol:displayId')
            sbol_displayid.text = component_displayid
            sbol_persistent = ET.SubElement(comp_xml, 'sbol:persistentIdentity')
            sbol_persistent.set('rdf:resource', compliant_component_uri)
            sbol_definition = ET.SubElement(comp_xml, 'sbol:definition')
            sbol_definition.set('rdf:resource', component_uri)
            try:
                comp_in_db = Feature.objects.get(uniquename=component_displayid)
            except Feature.DoesNotExist:
                comp_in_db = None

            sbol_access = ET.SubElement(comp_xml, 'sbol:access')
            if comp_in_db and comp_in_db.is_public:
                access_uri = 'http://sbols.org/v2#public'
            else:
                access_uri = 'http://sbols.org/v2#private'
            sbol_access.set('rdf:resource', access_uri)

#     if compliant_component_uris:
#         consts = ET.SubElement(comp_def_xml, "sbol:sequenceConstraint")
#         for index, subject_uri in enumerate(compliant_component_uris):
#             try:
#                 object_uri = compliant_component_uris[index + 1]
#             except IndexError:
#                 break
#             print subject_uri, object_uri
#             const = ET.SubElement(consts, "sbol:SequenceConstraint")
#             conts_uri = uri + '/r_' + str(index)
#             const.set('rdf:about', conts_uri)
#             ET.SubElement(const, 'sbol:object', {'rdf:resource': object_uri})
#             ET.SubElement(const, 'sbol:subject', {'rdf:resource': subject_uri})
#             ET.SubElement(const, 'sbol:restriction',
#                           {'rdf:resource': "http://sbols.org/v2#precedes"})
#     # get_role
#     try:
#         part_in_db = Feature.objects.get(uniquename=name)
#     except Feature.DoesNotExist:
#         part_in_db = None
#     if part_in_db:
#         type_ = part_in_db.type.name
#         so_uri = 'prueba'  # GOLDEN_TO_SO[type_]
#     else:
#         so_uri = 'generic'
#     sbol_role = ET.SubElement(comp_def_xml, 'sbol:role')
#     sbol_role.set('rdf:resource', so_uri)
    return comp_def_xml


def populate_sequence_element(seq_xml, seqrecord):
    seq_id = seqrecord.name + '_seq'
    uri = HOMESPACE + '/' + seq_id
    seq_xml.set('rdf:about', uri)
    displayid = ET.SubElement(seq_xml, 'sbol:displayId')
    displayid.text = seq_id
    elements = ET.SubElement(seq_xml, 'sbol:elements')
    elements.text = str(seqrecord.seq)
    encoding = ET.SubElement(seq_xml, 'sbol:encoding')
    encoding.set('rdf:resource', "www.chem.qmul.ac.uk/iubmb/misc/naseq.html")
    persisten_ident = ET.SubElement(seq_xml, 'sbol:persistentIdentity')
    persisten_ident.set('rdf:resource', uri)

    return seq_xml
