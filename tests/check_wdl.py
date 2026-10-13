"""Check Terra file creation and cloud-to-local command path handling."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import WDL

ROOT = Path(__file__).resolve().parents[1]
WRITERS = {"write_lines", "write_tsv", "write_map", "write_json", "write_object", "write_objects"}


def find_writers(node):
    found = []
    if isinstance(node, WDL.Expr.Apply) and node.function_name in WRITERS:
        found.append((node.function_name, node.pos.line))
    for child in node.children:
        found.extend(find_writers(child))
    return found


def main():
    document = WDL.load(str(ROOT / "workflows/pair_ld.wdl"))
    assert document.wdl_version == "1.0"
    assert not find_writers(document.workflow), "Workflow-scope file-writing function"
    task = document.tasks[0]
    assert not find_writers(task.command), "File creation in engine-side command expressions"
    # Simulate the engine's File localization step, then execute the rendered
    # command with an argument recorder. No genotype analysis or container runs.
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        recorder = work / "python"
        recorder.write_text('#!/usr/bin/env python3\nimport json, os, sys\nwith open(os.environ["ARG_LOG"], "w") as out:\n    json.dump(sys.argv[1:], out)\n')
        recorder.chmod(0o755)
        mapping = {}
        for present in (False, True):
            env = WDL.Env.Bindings()
            for name in ("pgen", "pvar", "psam", "variant_pairs", "sample_populations", "pgi"):
                cloud = f"gs://test-bucket/{name}"
                local = work / (name + " '$(touch INJECTED)")
                local.write_text("localized input\n")
                mapping[cloud] = str(local)
                env = env.bind(name, WDL.Value.File(cloud) if name != "pgi" or present else WDL.Value.Null())
            env = env.bind("cpu", WDL.Value.Int(2)).bind("plink_memory_mb", WDL.Value.Int(6000)).bind("allow_small_groups", WDL.Value.Boolean(present))
            localized = WDL.Value.rewrite_env_paths(env, lambda file: mapping[file.value])
            command = task.command.eval(localized, WDL.StdLib.Base("1.0")).value
            subprocess.run(["bash", "-c", command], check=True, cwd=work, capture_output=True, env={**os.environ, "PATH": str(work) + os.pathsep + os.environ["PATH"], "ARG_LOG": str(work / "arguments.json")})
            args = json.loads((work / "arguments.json").read_text())
            for flag, name in [("--pgen", "pgen"), ("--pvar", "pvar"), ("--psam", "psam"), ("--pairs", "variant_pairs"), ("--populations", "sample_populations")]:
                assert args[args.index(flag) + 1] == mapping[f"gs://test-bucket/{name}"]
            assert ("--pgi" in args) == present
            if present:
                assert args[args.index("--pgi") + 1] == mapping["gs://test-bucket/pgi"]
            assert ("--allow-small-groups" in args) == present
            assert not (work / "INJECTED").exists(), "Unsafe command quoting"
    print("PASS: WDL 1.0; no workflow-scope file writers; localized CLI paths and quoting")


if __name__ == "__main__":
    main()
