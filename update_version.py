#!/usr/bin/env python
"""
Helper functions for updating version
"""

import subprocess

import toml


def get_version():
    """
    Generate a version string based on the latest Git tag on main and the current branch.

    This function uses `git describe --tags main` to get the most recent tag on main and the number
    of commits since that tag. If the current branch is not `main`, it appends `-dev`
    to the version string. If no tag is found, it defaults to `0.1`.

    Returns:
        str: The generated version string.
    """
    try:
        # Describe latest tag on main with commit count
        describe = subprocess.check_output(["git", "describe", "--tags", "origin/main"]).strip().decode("utf-8")
        if "-" in describe:
            base, commits, _ = describe.split("-")
            major_minor = base.lstrip("v")  # Remove 'v' prefix if present
            version = f"{major_minor}.{commits}"
        else:
            version = describe.lstrip("v") + ".0"
    except subprocess.CalledProcessError:
        version = "0.1.0"

    # Get the current branch name
    branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]).strip().decode("utf-8")
    if branch != "main":
        version += "-dev"

    return version


def update_pyproject_version():
    """
    Update the version field in the `pyproject.toml` file based on the current Git state.

    This function generates a version string using `get_version()` and then updates the
    `version` field in the `pyproject.toml` file. The updated version is then saved back
    to the file.
    """
    version = get_version()

    # Load pyproject.toml
    with open("pyproject.toml", "r", encoding="utf-8") as file:
        pyproject_data = toml.load(file)

    # Update the version field
    pyproject_data["project"]["version"] = version

    # Write back the changes to pyproject.toml
    with open("pyproject.toml", "w", encoding="utf-8") as file:
        toml.dump(pyproject_data, file)

    print(f"Updated pyproject.toml version to {version}")


if __name__ == "__main__":
    update_pyproject_version()

    # # Get the current branch name
    # branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]).strip().decode("utf-8")
    # if branch == "main":
    #     # Round up version to the nearest decimal
    #     version_parts = version.split('.')
    #     if len(version_parts) >= 2:
    #         major = version_parts[0]
    #         minor = version_parts[1]
    #         try:
    #             new_minor = int(minor) + 1
    #             version = f"{major}.{new_minor}"
    #         except ValueError:
    #             # If minor is not an integer, keep the version as is
    #             pass
    #     else:
    #         # If version is like '0', increment it to '1'
    #         try:
    #             new_major = int(version_parts[0]) + 1
    #             version = f"{new_major}"
    #         except ValueError:
    #             # If major is not an integer, keep the version as is
    #             pass
    # else:
    #     # Append '-dev' for non-main branches
    #     version += "-dev"

    # return version
