#!/usr/bin/env python3
"""
Build Scripts for DEGAS Tools

Combine the source files in the `dev/` directory with `mtl_utils.py` to generate the script for the production environment
"""

import os
from pathlib import Path
import re
from datetime import datetime


class ScriptBuilder:
    """
    Construct the DEGAS tool script
    """

    def __init__(self, dev_dir, output_dir):
        self.dev_dir = Path(dev_dir)
        self.output_dir = Path(output_dir)
        self.utils_file = self.dev_dir / "mtl_utils.py"

        # read utils file
        with open(self.utils_file, "r", encoding="utf-8") as f:
            utils_content = f.read()

        # clean utils docstring
        self.utils_code = self._clean_utils_content(utils_content)

    def _clean_utils_content(self, content):
        """
        Clean up the contents of the utils file by removing the top-level docstring
        while retaining the import statements
        """

        lines = content.split("\n")
        result_lines = []
        in_docstring = False
        docstring_removed = False
        skip_empty = True  # Skip the initial blank lines.

        for line in lines:
            if skip_empty:
                if line.strip() == "":
                    continue
                skip_empty = False

            # skip the first docstring
            if not docstring_removed and '"""' in line:
                if not in_docstring:
                    in_docstring = True
                    continue
                else:
                    in_docstring = False
                    docstring_removed = True
                    continue

            if in_docstring and not docstring_removed:
                continue

            result_lines.append(line)

        # Remove the excess blank lines at the beginning
        while result_lines and result_lines[0].strip() == "":
            result_lines.pop(0)

        return "\n".join(result_lines)

    def _read_source_without_inline_marker(self, source_path):
        """
        Read the source file and remove the line marked with `INLINE_UTILS_HERE`
        
        NOTE: docstring is included in the result
        """

        with open(source_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        filtered_lines = [line for line in lines if "INLINE_UTILS_HERE" not in line]

        return "".join(filtered_lines)
    
    def build_script(self, source_file, output_name):
        """
        Construct one script file

        Args:
            source_file: source filename (e.g. 'BlankClassMTL_p3_src.py')
            output_name: output filename (e.g. 'BlankClassMTL_p3.py')
        """

        source_path = self.dev_dir / source_file
        output_path = self.output_dir / output_name

        if not source_path.exists():
            print(f"Warning: Source file not found: {source_path}")
            return False

        # Read the source file and remove the line marked with `INLINE_UTILS_HERE`
        source_content = self._read_source_without_inline_marker(source_path)

        # result
        parts = []

        # 1. Add an auto-generated docstring
        parts.append(self._generate_header(source_file))

        # 2. For the p1 file, inline the utils code no matter whether INLINE_UTILS_HERE is present.
        if "_p1_" in source_file:
            parts.append("\n")
            parts.append("# " + "=" * 70)
            parts.append("# Utility Functions (from mtl_utils.py)")
            parts.append("# " + "=" * 70)
            parts.append("\n")
            parts.append(self.utils_code)
            parts.append("\n\n")
            parts.append("# " + "=" * 70)
            parts.append("# Main Script")
            parts.append("# " + "=" * 70)
            parts.append("\n")

        # 3. Add the source code
        parts.append(source_content)

        output_content = "\n".join(parts)

        with open(output_path, "w", encoding="utf-8") as f:
            f.write(output_content)

        print(f"✓ Built: {output_name}")
        return True

    def _extract_docstring(self, filepath):
        """
        Retrieve the docstring from the source file (if any).
        
        This function is deprecated 
        """

        with open(filepath, "r", encoding="utf-8") as f:
            lines = f.readlines()

        docstring_lines = []
        in_docstring = False
        found_docstring = False

        for line in lines:
            if '"""' in line and not found_docstring:
                if not in_docstring:
                    in_docstring = True
                    docstring_lines.append(line)
                else:
                    docstring_lines.append(line)
                    found_docstring = True
                    break
            elif in_docstring:
                docstring_lines.append(line)

        if found_docstring:
            return "".join(docstring_lines)
        return None

    def _generate_header(self, source_file):
        """
        Generate the file header comment
        """
        
        return f"""
# {'=' * 70}
# Auto-generated from: {source_file}
# Build time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
# 
# DO NOT EDIT THIS FILE DIRECTLY
# Edit source files in the dev/ folder and rebuild using build_scripts.py
# {'=' * 70}

"""

    def build_all(self):
        """
        Construct all scripts
        """

        scripts = [
            # (srouce filename, output filename)
            ("ClassClassMTL_p1_src.py", "ClassClassMTL_p1.py"),
            ("ClassClassMTL_p3_src.py", "ClassClassMTL_p3.py"),
            ("ClassCoxMTL_p1_src.py", "ClassCoxMTL_p1.py"),
            ("ClassCoxMTL_p3_src.py", "ClassCoxMTL_p3.py"),
            ("BlankClassMTL_p1_src.py", "BlankClassMTL_p1.py"),
            ("BlankClassMTL_p3_src.py", "BlankClassMTL_p3.py"),
            ("BlankCoxMTL_p1_src.py", "BlankCoxMTL_p1.py"),
            ("BlankCoxMTL_p3_src.py", "BlankCoxMTL_p3.py"),
            ("ClassBlankMTL_p1_src.py", "ClassBlankMTL_p1.py"),
            ("ClassBlankMTL_p3_src.py", "ClassBlankMTL_p3.py")
        ]

        print("Building DEGAS Tools scripts...")
        print("=" * 70)
        print(f"Source directory: {self.dev_dir}")
        print(f"Output directory: {self.output_dir}")
        print(f"Utils file: {self.utils_file}")
        print("=" * 70)

        success_count = 0
        failed = []

        for source, output in scripts:
            if self.build_script(source, output):
                success_count += 1
            else:
                failed.append(source)

        print("=" * 70)
        print(f"Build complete: {success_count}/{len(scripts)} scripts generated")

        if failed:
            print("\nFailed to build:")
            for f in failed:
                print(f"  ✗ {f}")

        return success_count == len(scripts)


def main():
    script_dir = Path(__file__).parent
    dev_dir = script_dir
    output_dir = script_dir.parent

    print("\n" + "=" * 70)
    print("DEGAS Tools Build Script")
    print("=" * 70 + "\n")

    utils_file = dev_dir / "mtl_utils.py"
    if not utils_file.exists():
        print(f"ERROR: mtl_utils.py not found at {utils_file}")
        return 1

    builder = ScriptBuilder(dev_dir, output_dir)
    
    success = builder.build_all()

    print("\n" + "=" * 70)
    if success:
        print("SUCCESS: All scripts built successfully!")
    else:
        print("ERROR: Some scripts failed to build")
    print("=" * 70 + "\n")

    return 0 if success else 1


if __name__ == "__main__":
    main()
    # Diary:
    # BlankClassMTL test pass
    # BlankCoxMTL test pass
    # ClassClassMTL test pass