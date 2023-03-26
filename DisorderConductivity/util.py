import git

def get_git_hash():
  """Get git version hash of the head.

  Returns
  -------
  length-6 hexsha of the head, and 'unknown' if it fails.
  """
  try:
    repo = git.Repo(os.path.dirname(os.path.abspath(__file__)))
    return repo.git.rev_parse( repo.head.object.hexsha, short=6)
  except:
    return 'unknown'