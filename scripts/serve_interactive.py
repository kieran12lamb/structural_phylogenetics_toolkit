#!/usr/bin/env python3
"""Lightweight Interactive Tree Server & Pipeline Execution Bridge.

Serves the interactive tree viewer on http://localhost:8000 and provides
REST API endpoints to:
1. Proxy queries to the Viro3D API (bypassing browser CORS limits).
2. Trigger pipeline runs directly from the HTML interface.
3. Stream pipeline execution logs in real time to the in-browser console.
"""

import http.server
import json
import os
import subprocess
import sys
import threading
import urllib.parse
import urllib.request
import webbrowser
from pathlib import Path

PORT = 8000
REPO_DIR = Path(__file__).resolve().parent.parent

# Shared pipeline process state
pipeline_state = {
    "is_running": False,
    "process": None,
    "logs": [],
    "returncode": None,
    "command": ""
}
pipeline_lock = threading.Lock()


def run_pipeline_thread(cmd_args):
    global pipeline_state
    with pipeline_lock:
        pipeline_state["is_running"] = True
        pipeline_state["logs"] = [f"[System] Starting pipeline: {' '.join(cmd_args)}\n"]
        pipeline_state["returncode"] = None
        pipeline_state["command"] = " ".join(cmd_args)

    try:
        proc = subprocess.Popen(
            cmd_args,
            cwd=str(REPO_DIR),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        with pipeline_lock:
            pipeline_state["process"] = proc

        for line in proc.stdout:
            with pipeline_lock:
                pipeline_state["logs"].append(line)

        proc.wait()
        with pipeline_lock:
            pipeline_state["is_running"] = False
            pipeline_state["returncode"] = proc.returncode
            if proc.returncode == 0:
                pipeline_state["logs"].append("\n[System] Pipeline completed successfully!\n")
            else:
                pipeline_state["logs"].append(f"\n[System] Pipeline exited with status {proc.returncode}.\n")

    except Exception as e:
        with pipeline_lock:
            pipeline_state["is_running"] = False
            pipeline_state["logs"].append(f"\n[System Error] Failed to execute pipeline: {e}\n")


class InteractiveBridgeHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(REPO_DIR), **kwargs)

    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)

        if parsed.path in ("/", "/index.html"):
            self.path = "/interactive_tree.html"
            return super().do_GET()

        if parsed.path == "/api/status":
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()
            with pipeline_lock:
                data = {
                    "is_running": pipeline_state["is_running"],
                    "returncode": pipeline_state["returncode"],
                    "logs": pipeline_state["logs"][-60:],  # Return latest 60 lines
                    "command": pipeline_state["command"]
                }
            self.wfile.write(json.dumps(data).encode("utf-8"))
            return

        if parsed.path == "/api/viro3d_query":
            qs = urllib.parse.parse_qs(parsed.query)
            qualifier = qs.get("qualifier", ["glycoprotein"])[0]
            count = qs.get("count", ["10"])[0]

            try:
                v_url = f"https://viro3d.cvr.gla.ac.uk/api/proteins/protein_name/?qualifier={urllib.parse.quote(qualifier)}&page_size={count}&page_num=1"
                req = urllib.request.Request(v_url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=15) as r:
                    v_data = json.loads(r.read().decode("utf-8"))

                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(json.dumps(v_data).encode("utf-8"))
            except Exception as e:
                self.send_response(500)
                self.send_header("Content-Type", "application/json")
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(json.dumps({"error": str(e)}).encode("utf-8"))
            return

        return super().do_GET()

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)

        if parsed.path == "/api/run_pipeline":
            content_len = int(self.headers.get("Content-Length", 0))
            body = self.rfile.read(content_len).decode("utf-8")
            try:
                payload = json.loads(body)
            except Exception:
                payload = {}

            with pipeline_lock:
                if pipeline_state["is_running"]:
                    self.send_response(409)
                    self.send_header("Content-Type", "application/json")
                    self.send_header("Access-Control-Allow-Origin", "*")
                    self.end_headers()
                    self.wfile.write(json.dumps({"error": "A pipeline job is already running."}).encode("utf-8"))
                    return

            # Construct command args
            cmd = [sys.executable, str(REPO_DIR / "scripts/viral_phylogenetics.py"), "pipeline"]
            
            source = payload.get("source", "viro3d")
            if source == "local":
                in_folder = payload.get("input_folder", "structures")
                cmd.extend(["-i", in_folder])
                if payload.get("metadata"):
                    cmd.extend(["-meta", payload["metadata"]])
            else:
                qualifier = payload.get("qualifier", "glycoprotein")
                count = str(payload.get("count", 50))
                cmd.extend(["-q", qualifier, "-c", count])

            aligner = payload.get("aligner", "foldmason")
            if aligner == "mafft":
                cmd.extend(["--aligner", "mafft"])

            tree_type = payload.get("tree_type", "both")
            cmd.extend(["--tree-type", tree_type])

            matrix = payload.get("matrix", "alphafold")
            cmd.extend(["--matrix", matrix])

            out_dir = payload.get("output_dir", "viral_custom_workflow")
            cmd.extend(["-o", out_dir])

            threads = str(payload.get("threads", "AUTO"))
            cmd.extend(["-t", threads])

            if payload.get("fast", False):
                cmd.append("--fast")
            elif payload.get("bootstrap", True):
                cmd.extend(["-b", "1000", "--alrt", "1000"])
            else:
                cmd.extend(["-b", "0", "--alrt", "0"])

            if payload.get("embed", False):
                cmd.append("--embed")
                if payload.get("embed_model"):
                    cmd.extend(["--embed-model", str(payload["embed_model"])])
                if payload.get("embed_clustering"):
                    cmd.extend(["--embed-clustering", str(payload["embed_clustering"])])

            # Start thread
            t = threading.Thread(target=run_pipeline_thread, args=(cmd,), daemon=True)
            t.start()

            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()
            self.wfile.write(json.dumps({"status": "started", "command": " ".join(cmd)}).encode("utf-8"))
            return

        self.send_response(404)
        self.end_headers()


def main():
    server = http.server.ThreadingHTTPServer(("0.0.0.0", PORT), InteractiveBridgeHandler)
    url = f"http://localhost:{PORT}"
    print("=" * 70)
    print(f"🚀 Viral Structural Phylogenetics Interactive Server Bridge")
    print(f"📡 Serving: {url}")
    print(f"📂 Workspace: {REPO_DIR}")
    print(f"✨ Features:")
    print(f"   - Interactive SVG & WebGL Phylogenetic Suite with 3D Structure Viewer")
    print(f"   - Live Viro3D API query proxy")
    print(f"   - Direct in-browser pipeline execution & real-time log console")
    print("=" * 70)

    try:
        webbrowser.open(url)
    except Exception:
        pass

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[Server] Shutting down.")
        server.server_close()


if __name__ == "__main__":
    main()
