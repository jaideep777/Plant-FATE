
**
TO BE UPDATED
**

### Updating external libraries (Phydro, libpspm, and Flare)

Do not edit the code in external/ folder unless you are familiar with git subtrees and know how to contribute it back to the parent repos. 

#### Adding a new external lib from a specific ref

git remote add -f `newlib` https://github.com/xxx/newlib
git subtree add --prefix external/`newlib` `newlib` `ref` --squash

#### updating an external library from a specific ref

git subtree pull --prefix external/`lib` `lib` `ref` --squash

#### *Danger Zone:* Updating the parent library with changes made locally

git subtree push --prefix external/`lib` `lib` `ref`
