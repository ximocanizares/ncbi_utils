import json
from pathlib import Path
import pickle
import hashlib
import xml.etree.ElementTree as ET

import requests

NCBI_EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_SEARCH_BASE_URL = NCBI_EUTILS_BASE_URL + "esearch.fcgi?"
NCBI_FETCH_BASE_URL = NCBI_EUTILS_BASE_URL + "efetch.fcgi?"


def _hash(value):
    try:
        hash(value)
    except TypeError:
        raise ValueError(
            "All arguments should be hasheable, but this one failed: " + str(value)
        )
    pickled = pickle.dumps(value)
    return hashlib.md5(pickled).hexdigest()


def cache_call(funct, cache_dir: Path, args=None, kwargs=None):
    if args is None:
        args = tuple()
    if kwargs is None:
        kwargs = {}

    hashes = [_hash(arg) for arg in args]
    hashes.extend((_hash(kwargs[arg]) for arg in sorted(kwargs.keys())))

    hash_ = _hash(tuple(hashes))
    cache_dir.mkdir(exist_ok=True)
    cache_dir = cache_dir / funct.__name__
    cache_dir.mkdir(exist_ok=True)
    cache_path = cache_dir / str(hash_)
    if cache_path.exists():
        with cache_path.open("rb") as fhand:
            result = pickle.load(fhand)
    else:
        result = funct(*args, **kwargs)
        with cache_path.open("wb") as fhand:
            pickle.dump(result, fhand)
    return result


def search_id_for_experiment_acc(acc: str) -> str:
    url = NCBI_SEARCH_BASE_URL + f"db=sra&term={acc}[Accession]&retmode=json&retmax=1"
    return _search_id_with(url, acc=acc, db="sra")


def _search_id_with(url, acc, db) -> str:
    response = requests.get(url)
    assert response.status_code == 200
    content = response.content
    search_result = json.loads(content)
    if not search_result["esearchresult"]:
        raise ValueError(f"acc {acc} not found in the {db} database")

    if "idlist" not in search_result["esearchresult"]:
        raise RuntimeError(f"Missing idlist for url: {url}")

    if not search_result["esearchresult"]["idlist"]:
        raise RuntimeError(f"No ids in idlist for acc: {acc}")
    elif len(search_result["esearchresult"]["idlist"]) != 1:
        raise RuntimeError(f"We expected just 1 id in idlist for acc: {acc}")

    id_ = search_result["esearchresult"]["idlist"][0]
    if not id_.isdigit():
        raise RuntimeError(
            f"We expected an all digit id for acc {acc}, but we got: {id_}"
        )

    return id_


def search_id_for_biosample_acc(acc: str) -> str:
    db = "biosample"
    url = NCBI_SEARCH_BASE_URL + f"db={db}&term={acc}[Accession]&retmode=json&retmax=1"
    return _search_id_with(url, acc, db=db)


def search_id_for_bioproject_acc(acc: str) -> str:
    url = (
        NCBI_SEARCH_BASE_URL
        + f"db=bioproject&term={acc}[Accession]&retmode=json&retmax=1"
    )
    return _search_id_with(url, acc, db="bioproject")


def fetch_bioproject_acc_for_experiment(id: str):
    if not id.isdigit():
        raise RuntimeError(f"We expected an all digit experiment id, but we got: {id}")

    url = NCBI_FETCH_BASE_URL + f"db=sra&id={id}"
    response = requests.get(url)
    assert response.status_code == 200
    xml = response.content
    experiment_set = ET.fromstring(xml)
    external_id = (
        experiment_set.find("EXPERIMENT_PACKAGE")
        .find("STUDY")
        .find("IDENTIFIERS")
        .find("EXTERNAL_ID")
    )
    if external_id is not None:
        return external_id.text
    else:
        raise RuntimeError(f"No bioproject acc found for experiment: {id}")


def fetch_bioproject_info(id: str):
    if not id.isdigit():
        raise RuntimeError(f"We expected an all digit experiment id, but we got: {id}")
    url = NCBI_FETCH_BASE_URL + f"db=bioproject&id={id}"

    response = requests.get(url)

    if response.status_code != 200:
        raise RuntimeError(f"There was an error requesting the url: {url}")

    xml = response.content
    record_set = ET.fromstring(xml)
    project = record_set.find("DocumentSummary").find("Project")
    info = {}
    info["acc"] = project.find("ProjectID").find("ArchiveID").attrib["accession"]
    info["title"] = project.find("ProjectDescr").find("Title").text
    info["description"] = project.find("ProjectDescr").find("Description").text
    try:
        target = (
            project.find("ProjectType").find("ProjectTypeSubmission").find("Target")
        )
    except AttributeError:
        target = None
    if target:
        info["type"] = {
            "capture": target.attrib["capture"],
            "material": target.attrib["material"],
        }
    else:
        info["type"] = {
            "capture": "",
            "material": "",
        }

    try:
        info["submission_last_update"] = (
            record_set.find("DocumentSummary").find("Submission").attrib["last_update"]
        )
    except KeyError:
        pass
    return info


def ask_ncbi_for_biosample_ids_in_bioproject(bioproject_id: str) -> list[str]:
    if bioproject_id.lower().startswith("prj"):
        raise ValueError(
            f"Use the numeric id, not the PRJXXXX accession: {bioproject_id}"
        )
    try:
        int(bioproject_id)
    except ValueError:
        raise ValueError(
            f"The id should be an str, but all numbers (e.g. 1025377), but it was: {bioproject_id}"
        )

    url = f"{NCBI_EUTILS_BASE_URL}elink.fcgi?dbfrom=bioproject&db=biosample&id={bioproject_id}&retmode=json"
    response = requests.get(url)

    if response.status_code != 200:
        raise RuntimeError(f"Thre was a problem getting the URL: {url}")

    jsons = response.content

    search_result = json.loads(jsons)

    biosample_ids = set()
    for linkset in search_result["linksets"]:
        for linksetdb in linkset["linksetdbs"]:
            biosample_ids.update(linksetdb["links"])
    return sorted(biosample_ids)


def fetch_biosample_info(biosample_id: str):
    if biosample_id.lower().startswith("prj"):
        raise ValueError(
            f"Use the numeric id, not the SAMNXXXX accession: {biosample_id}"
        )
    try:
        int(biosample_id)
    except ValueError:
        raise ValueError(
            f"The id should be an str, but all numbers (e.g. 1025377), but it was: {biosample_id}"
        )

    url = f"{NCBI_EUTILS_BASE_URL}efetch.fcgi?db=biosample&id={biosample_id}"
    response = requests.get(url)

    if response.status_code != 200:
        raise RuntimeError(f"There was a problem getting the URL: {url}")

    xml = response.content

    biosample_set = ET.fromstring(xml)
    biosample_xml = biosample_set.find("BioSample")
    biosample = {}
    biosample["biosampledb_id"] = biosample_xml.attrib["id"]
    biosample["biosampledb_accession"] = biosample_xml.attrib["accession"]
    biosample["publication_date"] = biosample_xml.attrib["publication_date"]
    for id in biosample_xml.find("Ids").findall("Id"):
        try:
            id.attrib["db"]
        except KeyError:
            continue
        if id.attrib["db"] == "SRA":
            biosample["sra_accession"] = id.text

    description = biosample_xml.find("Description")
    biosample["title"] = description.find("Title").text
    organism = description.find("Organism")
    biosample["organism_id"] = organism.attrib["taxonomy_id"]
    biosample["organism_name"] = organism.attrib["taxonomy_name"]

    attributes = {}
    for attribute_xml in biosample_xml.find("Attributes"):
        attributes[attribute_xml.attrib["attribute_name"]] = attribute_xml.text
    biosample["attributes"] = attributes

    return biosample


def search_experiments_in_sra_with_biosample_accession(biosample_acc, cache_dir=None):
    url = f"{NCBI_EUTILS_BASE_URL}esearch.fcgi?db=sra&term={biosample_acc}[BioSample]&retmode=json"
    response = requests.get(url)
    assert response.status_code == 200
    jsons = response.content
    search_result = json.loads(jsons)
    ids = search_result["esearchresult"]["idlist"]

    return ids


def fetch_experiment_info(experiment_id: str):
    if experiment_id.lower().startswith("s"):
        raise ValueError(
            f"Use the numeric id, not the SRXXXXX accession: {experiment_id}"
        )
    try:
        int(experiment_id)
    except ValueError:
        raise ValueError(
            f"The id should be an str, but all numbers (e.g. 1025377), but it was: {experiment_id}"
        )

    url = f"{NCBI_EUTILS_BASE_URL}efetch.fcgi?db=sra&id={experiment_id}"
    response = requests.get(url)
    xml = response.content
    experiment_package_set = ET.fromstring(xml)
    experiment_packages = experiment_package_set.findall("EXPERIMENT_PACKAGE")
    if not experiment_packages:
        raise RuntimeError(f"No experiment package found: {url}")
    elif len(experiment_packages) != 1:
        raise RuntimeError(f"We expected only one experiment package: {url}")
    experiment_package = experiment_packages[0]
    experiments = experiment_package.findall("EXPERIMENT")
    if len(experiments) != 1:
        raise RuntimeError(f"We expected only one experiment: {url}")
    experiment = experiments[0]

    info = {"acc": experiment.find("IDENTIFIERS").find("PRIMARY_ID").text}
    info["study_acc"] = experiment.find("STUDY_REF").attrib["accession"]
    info["Bioproject_acc"] = experiment.find("STUDY_REF").find("IDENTIFIERS").find("EXTERNAL_ID").text
    
    info["design"] = {
        "description": experiment.find("DESIGN").find("DESIGN_DESCRIPTION").text,
        "sample": experiment.find("DESIGN")
        .find("SAMPLE_DESCRIPTOR")
        .attrib["accession"],
        "library": {
            "strategy": experiment.find("DESIGN")
            .find("LIBRARY_DESCRIPTOR")
            .find("LIBRARY_STRATEGY")
            .text,
            "source": experiment.find("DESIGN")
            .find("LIBRARY_DESCRIPTOR")
            .find("LIBRARY_SOURCE")
            .text,
            "selection": experiment.find("DESIGN")
            .find("LIBRARY_DESCRIPTOR")
            .find("LIBRARY_SELECTION")
            .text,
            "layout": experiment.find("DESIGN")
            .find("LIBRARY_DESCRIPTOR")
            .find("LIBRARY_LAYOUT")
            .text,
        },
    }
    platform = experiment.find("PLATFORM")
    if platform.find("ILLUMINA"):
        info["platform"] = {
            "platform": "ilumina",
            "instrument_model": platform.find("ILLUMINA").find("INSTRUMENT_MODEL").text,
        }
    elif platform.find("DNBSEQ"):
        info["platform"] = {
            "platform": "dnbseq",
            "instrument_model": platform.find("DNBSEQ").find("INSTRUMENT_MODEL").text,
        }
    elif platform.find("OXFORD_NANOPORE"):
        info["platform"] = {
            "platform": "nanopore",
            "instrument_model": platform.find("OXFORD_NANOPORE")
            .find("INSTRUMENT_MODEL")
            .text,
        }
    elif platform.find("PACBIO_SMRT"):
        info["platform"] = {
            "platform": "pacbio",
            "instrument_model": platform.find("PACBIO_SMRT")
            .find("INSTRUMENT_MODEL")
            .text,
        }
    elif platform.find("BGISEQ"):
        info["platform"] = {
            "platform": "bgiseq",
            "instrument_model": platform.find("BGISEQ").find("INSTRUMENT_MODEL").text,
        }
    elif platform.find("LS454"):
        info["platform"] = {
            "platform": "454",
            "instrument_model": platform.find("LS454").find("INSTRUMENT_MODEL").text,
        }
    elif platform.find("ABI_SOLID"):
        info["platform"] = {
            "platform": "solid",
            "instrument_model": platform.find("ABI_SOLID").find("INSTRUMENT_MODEL").text,
        }
    elif platform.find("ELEMENT"):
        info["platform"] = {
            "platform": "element",
            "instrument_model": platform.find("ELEMENT").find("INSTRUMENT_MODEL").text,
        }
    else:
        raise RuntimeError(f"Unknown platform for: {url}")

    info["organization"] = experiment_package.find("Organization").find("Name").text
    try:
        info["study"] = (
            experiment_package.find("STUDY")
            .find("DESCRIPTOR")
            .find("STUDY_DESCRIPTION")
            .text
        )
    except AttributeError:
        info["study"] = (
            experiment_package.find("STUDY").find("DESCRIPTOR").find("STUDY_TITLE").text
        )

    runs = []
    run_set = experiment_package.find("RUN_SET")
    for run in run_set.findall("RUN"):
        run_info = {}
        run_info["accession"] = run.attrib["accession"]
        run_info["date"] = run.attrib["published"]
        run_info["is_public"] = run.attrib["is_public"]

        files = []
        try:
            sra_files = run.find("SRAFiles").findall("SRAFile")
        except AttributeError:
            sra_files = []
        for file in sra_files:
            file_info = {}
            file_info["name"] = file.attrib["filename"]
            file_info["md5"] = file.attrib["md5"]
            file_info["s3_url"] = file.find("Alternatives").attrib["url"]
            files.append(file_info)
        run_info["files"] = files

        runs.append(run_info)
    info["runs"] = runs

    return info


if __name__ == "__main__":
    # id_ = search_id_for_experiment_acc("SRX27341610")
    cache_dir = Path("__file__").absolute().parent / "cache"
    cache_dir.mkdir(exist_ok=True)

    experiments_to_ignore = ["SRX118405"]
    bioprojects = [
        "PRJEB29506",
        "PRJEB30666",
        "PRJEB31935",
        "PRJEB43158",
        "PRJEB43865",
        "PRJEB4395",
        "PRJEB44956",
        "PRJEB4585",
        "PRJEB59321",
        "PRJEB6302",
        "PRJEB63089",
        "PRJEB80388",
        "PRJEB9878",
        "PRJNA1006521",
        "PRJNA1007367",
        "PRJNA1066603",
        "PRJNA1078922",
        "PRJNA1090062",
        "PRJNA1100794",
        "PRJNA1108418",
        "PRJNA1108922",
        "PRJNA1128095",
        "PRJNA1162684",
        "PRJNA1193526",
        "PRJNA1210015",
        "PRJNA212546",
        "PRJNA222545",
        "PRJNA239450",
        "PRJNA240684",
        "PRJNA240686",
        "PRJNA240687",
        "PRJNA240689",
        "PRJNA240691",
        "PRJNA240692",
        "PRJNA240693",
        "PRJNA240694",
        "PRJNA240695",
        "PRJNA240696",
        "PRJNA240698",
        "PRJNA240699",
        "PRJNA240701",
        "PRJNA240702",
        "PRJNA240703",
        "PRJNA240704",
        "PRJNA240707",
        "PRJNA240708",
        "PRJNA240709",
        "PRJNA240712",
        "PRJNA240763",
        "PRJNA240764",
        "PRJNA240776",
        "PRJNA240777",
        "PRJNA240779",
        "PRJNA240780",
        "PRJNA240781",
        "PRJNA240782",
        "PRJNA240783",
        "PRJNA240785",
        "PRJNA240786",
        "PRJNA240787",
        "PRJNA240788",
        "PRJNA240790",
        "PRJNA240791",
        "PRJNA240792",
        "PRJNA240793",
        "PRJNA240794",
        "PRJNA240796",
        "PRJNA240797",
        "PRJNA240801",
        "PRJNA240803",
        "PRJNA240804",
        "PRJNA240805",
        "PRJNA240806",
        "PRJNA240807",
        "PRJNA240809",
        "PRJNA240810",
        "PRJNA240811",
        "PRJNA240813",
        "PRJNA240814",
        "PRJNA240816",
        "PRJNA240817",
        "PRJNA240818",
        "PRJNA240821",
        "PRJNA240823",
        "PRJNA240824",
        "PRJNA240825",
        "PRJNA240826",
        "PRJNA240827",
        "PRJNA240828",
        "PRJNA240829",
        "PRJNA240830",
        "PRJNA240833",
        "PRJNA240834",
        "PRJNA240835",
        "PRJNA240836",
        "PRJNA240838",
        "PRJNA240839",
        "PRJNA240840",
        "PRJNA240841",
        "PRJNA240842",
        "PRJNA240843",
        "PRJNA240844",
        "PRJNA240845",
        "PRJNA240846",
        "PRJNA240849",
        "PRJNA240851",
        "PRJNA240852",
        "PRJNA240853",
        "PRJNA240854",
        "PRJNA240855",
        "PRJNA240856",
        "PRJNA240857",
        "PRJNA240859",
        "PRJNA240860",
        "PRJNA240861",
        "PRJNA240862",
        "PRJNA240865",
        "PRJNA240866",
        "PRJNA240867",
        "PRJNA240868",
        "PRJNA240870",
        "PRJNA240917",
        "PRJNA240918",
        "PRJNA240919",
        "PRJNA240920",
        "PRJNA240921",
        "PRJNA240922",
        "PRJNA240923",
        "PRJNA240924",
        "PRJNA240925",
        "PRJNA240927",
        "PRJNA240928",
        "PRJNA240929",
        "PRJNA240931",
        "PRJNA240932",
        "PRJNA240933",
        "PRJNA240934",
        "PRJNA240936",
        "PRJNA240937",
        "PRJNA240938",
        "PRJNA240939",
        "PRJNA240940",
        "PRJNA240941",
        "PRJNA240942",
        "PRJNA240943",
        "PRJNA240944",
        "PRJNA240946",
        "PRJNA259308",
        "PRJNA295106",
        "PRJNA312569",
        "PRJNA313013",
        "PRJNA325096",
        "PRJNA353136",
        "PRJNA353161",
        "PRJNA356494",
        "PRJNA376115",
        "PRJNA378916",
        "PRJNA387314",
        "PRJNA391809",
        "PRJNA423252",
        "PRJNA429495",
        "PRJNA438044",
        "PRJNA447965",
        "PRJNA454805",
        "PRJNA494243",
        "PRJNA494932",
        "PRJNA505969",
        "PRJNA509653",
        "PRJNA527863",
        "PRJNA540313",
        "PRJNA540748",
        "PRJNA550259",
        "PRJNA557253",
        "PRJNA566320",
        "PRJNA574073",
        "PRJNA592130",
        "PRJNA615189",
        "PRJNA616074",
        "PRJNA629854",
        "PRJNA637170",
        "PRJNA643909",
        "PRJNA666021",
        "PRJNA668905",
        "PRJNA672142",
        "PRJNA673252",
        "PRJNA678238",
        "PRJNA692070",
        "PRJNA695390",
        "PRJNA697846",
        "PRJNA702633",
        "PRJNA703300",
        "PRJNA704807",
        "PRJNA708163",
        "PRJNA713664",
        "PRJNA725647",
        "PRJNA750735",
        "PRJNA751989",
        "PRJNA776037",
        "PRJNA776853",
        "PRJNA778961",
        "PRJNA778983",
        "PRJNA788715",
        "PRJNA790656",
        "PRJNA796560",
        "PRJNA802260",
        "PRJNA846963",
        "PRJNA863724",
        "PRJNA870056",
        "PRJNA882342",
        "PRJNA884656",
        "PRJNA902481",
        "PRJNA907318",
        "PRJNA931572",
        "PRJNA935397",
        "PRJNA936910",
        "PRJNA953768",
        "PRJNA988534",
        "PRJEB5235",
        "SRP010718",
    ]
    biosample_ids_for_bioprojects = {"SRP010718": ["780126"]}
    experiments = {}
    for idx, bioproject_acc in enumerate(bioprojects):
        if bioproject_acc in biosample_ids_for_bioprojects:
            biosample_ids = biosample_ids_for_bioprojects[bioproject_acc]
        else:
            project_id = cache_call(
                search_id_for_bioproject_acc,
                args=(bioproject_acc,),
                cache_dir=cache_dir,
            )
            project_info = cache_call(
                fetch_bioproject_info, args=(project_id,), cache_dir=cache_dir
            )
            biosample_ids = cache_call(
                ask_ncbi_for_biosample_ids_in_bioproject,
                args=(project_id,),
                cache_dir=cache_dir,
            )

        for biosample_id in biosample_ids:
            biosample = cache_call(
                fetch_biosample_info, args=(biosample_id,), cache_dir=cache_dir
            )
            experiment_ids = cache_call(
                search_experiments_in_sra_with_biosample_accession,
                args=(biosample["biosampledb_accession"],),
                cache_dir=cache_dir,
            )
            for experiment_id in experiment_ids:
                experiment = cache_call(
                    fetch_experiment_info, args=(experiment_id,), cache_dir=cache_dir
                )
                date = [run["date"] for run in experiment["runs"]][0]
                is_public = [run["date"] for run in experiment["runs"]][0]

                for run in experiment["runs"]:
                    print(run["accession"])

                experiments[experiment["acc"]] = {
                    "accession": experiment["acc"],
                    "bioproject": bioproject_acc,
                    "bioproject_title": project_info["title"],
                    "bioproject_description": project_info["description"],
                    "biosample_title": biosample["title"],
                    "biosample_organism": biosample["organism_name"],
                    "biosample_cultivar": biosample.get("attributes", {}).get(
                        "cultivar", ""
                    ),
                    "library_strategy": experiment["design"]["library"]["strategy"],
                    "library_source": experiment["design"]["library"]["source"],
                    "library_selection": experiment["design"]["library"]["selection"],
                    "platform": experiment["platform"]["instrument_model"],
                    "organization": experiment["organization"],
                    "study": experiment["study"],
                    "date": date,
                    "is_public": is_public,
                }

    import pandas

    experiments = pandas.DataFrame(experiments.values())
    print(experiments)
    experiments.to_excel("../../experiments.xlsx", index=False)
