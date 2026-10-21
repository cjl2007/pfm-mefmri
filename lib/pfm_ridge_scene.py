#!/usr/bin/env python3
"""Adapt the crawling-seed flat/inflated scene to a static RidgePFM map."""
from __future__ import annotations
import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from crawling_seed_fc_movie import rewrite_surface_references


def build_scene(template, cifti, surface_dir, flat_dir, subject, output, map_index=0, border_paths=(), compact=False):
    text=rewrite_surface_references(Path(template).read_text(),str(Path(surface_dir).resolve()),subject,str(Path(flat_dir).resolve()),subject)
    root=ET.fromstring(text)
    scenes=root.findall('Scene')
    for scene in scenes[1:]:root.remove(scene)
    infos=root.find('SceneInfoDirectory')
    if infos is not None:
        for info in infos.findall('SceneInfo')[1:]:infos.remove(info)
    cifti=Path(cifti).resolve();output=Path(output).resolve()
    # Remove unavailable custom palettes and saved connectivity-only state.
    changed_arrays=[]
    for parent in root.iter():
        for child in list(parent):
            if child.tag=='Object' and child.get('Class')=='ConnectivityDataLoaded':parent.remove(child)
            elif child.tag=='Element' and any(n.get('Class') == 'CaretDataFile' and n.get('Name','').endswith('.palette') for n in child.findall('Object')):
                parent.remove(child);changed_arrays.append(parent)
    for arr in changed_arrays:
        elements=arr.findall('Element')
        if elements:
            arr.set('Length',str(len(elements)))
            for i,e in enumerate(elements):e.set('Index',str(i))
    for node in root.iter('Object'):
        text=node.text or ''
        if '.dconn.nii' in text and node.get('Type') in ('pathName','string'):
            node.text=str(cifti) if node.get('Type')=='pathName' else cifti.name
        if node.get('Name','').endswith('.dconn.nii'):node.set('Name',cifti.name)
        if node.get('Name') in ('m_showSurfaceIdentificationSymbols','m_showVolumeIdentificationSymbols'):node.text='false'
        if node.get('Name')=='selectedMapName':node.text=''
        if node.get('Name')=='selectedMapIndex':node.text=str(map_index)
        if node.get('Type')=='pathName' and text.endswith('.scene'):node.text=str(output)
        if node.get('Name')=='m_simpleText' and ('CONNECTIVITY' in text or 'VERTEX' in text):node.text=''
    # File type in the embedded specification must match the selected overlay.
    for node in root.iter():
        if node.text and 'CONNECTIVITY_DENSE' in node.text:
            node.text=node.text.replace('CONNECTIVITY_DENSE','CONNECTIVITY_DENSE_LABEL' if cifti.name.endswith('.dlabel.nii') else 'CONNECTIVITY_DENSE_SCALAR')
    if compact:
        for node in root.iter('Object'):
            if node.get('Name') in ('m_customTileTabsConfigurationLatest','m_customTileTabsConfiguration','m_sceneTileTabsConfiguration'):
                node.text=(node.text or '').replace('WeightStretch="0.70"','WeightStretch="0.45"').replace('>1 0.7<','>1 0.45<')
    if border_paths:
        if len(border_paths)!=2:raise ValueError('Supply left and right border files')
        caret=root.find('.//ObjectArray[@Name="allCaretDataFiles_V2"]')
        spec=root.find('.//ObjectArray[@Name="dataFilesArray"]')
        for path,structure in zip(border_paths,['CORTEX_LEFT','CORTEX_RIGHT']):
            path=Path(path).resolve()
            if not path.is_file():raise FileNotFoundError(path)
            for array,cls,name,fields in [
                (caret,'CaretDataFile',path.name,[('pathName','dataFileName_V2',str(path))]),
                (spec,'SpecFileDataFile','specFileDataFile',[
                    ('enumeratedType','dataFileType','BORDER'),('enumeratedType','structure',structure),
                    ('pathName','fileName',str(path)),('boolean','selected','true')])]:
                element=ET.SubElement(array,'Element',Index=str(len(array)))
                obj=ET.SubElement(element,'Object',Type='class',Class=cls,Name=name,Version='1')
                for typ,key,value in fields:ET.SubElement(obj,'Object',Type=typ,Name=key).text=value
                array.set('Length',str(len(array)))
        display=root.find('.//Object[@Class="DisplayPropertiesBorders"]')
        for node in display:
            name=node.get('Name','')
            value=('true' if name.startswith('m_displayStatusIn') else
                   'FEATURE_COLORING_TYPE_STANDARD_COLOR' if name.startswith('m_coloringTypeIn') else
                   'DRAW_AS_LINES' if name.startswith('m_drawingTypeIn') else
                   'BLACK' if name.startswith('m_standardColorTypeIn') else None)
            if value is not None:
                for entry in node:entry.text=value
    output.parent.mkdir(parents=True,exist_ok=True)
    ET.ElementTree(root).write(output,encoding='UTF-8',xml_declaration=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['template','cifti','surface-dir','flat-dir','subject','output']:p.add_argument('--'+key,required=True)
    p.add_argument('--map-index',type=int,default=0);a=p.parse_args()
    build_scene(a.template,a.cifti,a.surface_dir,a.flat_dir,a.subject,a.output,a.map_index)
if __name__=='__main__':main()
